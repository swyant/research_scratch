
struct PACE_Ensemble
    yace_fnames::Array{String}
    lmps::Array{LMP}
    energy_thresh::Float64
    force_thresh::Float64 
    extra_data::Bool
end

function PACE_Ensemble(yace_fnames::Array{String}, elem_strs::Array{String}; energy_thresh=1.0,force_thresh=1.0,extra_data=false)
    lmp_list = []
    for fname in yace_fnames
        push!(lmp_list,initialize_committee_member(fname, elem_strs))
    end
    pace_ens = PACE_Ensemble(yace_fnames,lmp_list,energy_thresh,force_thresh,extra_data)

    pace_ens
end

function close!(pace_cmte::PACE_Ensemble)
    for c_lmps in pace_cmte.lmps
        LAMMPS.API.lammps_close(c_lmps)
    end
end

function initialize_committee_member(ace_fname, atom_type_list)
    num_types = length(atom_type_list)
    atomtype_str = ""
    for elem_str in atom_type_list
        atomtype_str = atomtype_str * " " * elem_str
    end
    print(atomtype_str)

    lmp = LMP(["-screen", "none"])
    command(lmp, "log none")

    command(lmp, "units          metal")
    command(lmp, "boundary       p p p")
    command(lmp, "atom_style     atomic")
    command(lmp, "neigh_modify   delay 0 every 1 check no") # probably not necessary for a single calc

    command(lmp, "region main block 0.0 1.0 0.0 1.0 0.0 1.0")
    command(lmp, "create_box $(num_types) main")

    command(lmp, "pair_style  pace")
    command(lmp, """pair_coeff     * * $(ace_fname) $(atomtype_str)""")

    command(lmp, "thermo       1")
    command(lmp, "thermo_style custom step temp pe ke etotal press")
    command(lmp,"compute pe all pe")
    return lmp
end

function initialize_uq_dict(ace_cmte::PACE_Ensemble)
    uq_dict = Dict("mean_fcomp_stdevs" => [],
                   "energy_stdevs" => [],
                   )
    
    if ace_cmte.extra_data 
        uq_dict["each_fcomp_mean"] = []
        uq_dict["each_fcomp_stdevs"] = []
        uq_dict["mean_pe"] = []
        uq_dict["cmte_pe_preds"] = []
    end

    uq_dict
end

function compute_force_comp_data(cmte_dict)
    f_stdevs = [] # 3N standard deviations 
    f_means  = []
    for f_idx in eachindex(cmte_dict[1]["forces"])
        tmp_fcomp_arr = [cmte_dict[cmte_key]["forces"][f_idx] for cmte_key in keys(cmte_dict)]
        f_comp_std = Statistics.std(tmp_fcomp_arr) 
        f_comp_mean = Statistics.mean(tmp_fcomp_arr)
        push!(f_stdevs, f_comp_std)
        push!(f_means, f_comp_mean)
    end
    avg_fcomp_stdev = mean(f_stdevs)
    avg_fcomp_stdev, transpose(reshape(f_stdevs,3,:)), transpose(reshape(f_means,3,:))
end

function compute_energy_mean_and_stdev(cmte_dict)
    e_stdev = Statistics.std([ cmte_dict[cmte_key]["pe"] for cmte_key in keys(cmte_dict)])
    e_mean = Statistics.mean([ cmte_dict[cmte_key]["pe"] for cmte_key in keys(cmte_dict)])
    e_mean, e_stdev
end


function update_uq!(ace_cmte::PACE_Ensemble, main_cfg_dict,uq_dict, uncertain_configs)
    cmte_data= Dict()
    cmte_lmps = ace_cmte.lmps
    for i in 1:length(cmte_lmps)
        cfg_dict = ace_singlepoint(main_cfg_dict["box_bounds"],main_cfg_dict["types"],main_cfg_dict["positions"],main_cfg_dict["masses"],0,cmte_lmps[i])
        cmte_data[i] = cfg_dict
    end
    avg_fcomp_stdev, fcomp_stdevs, fcomp_means = compute_force_comp_data(cmte_data)
    energy_mean, energy_stdev = compute_energy_mean_and_stdev(cmte_data)
    raw_energies = [cmte_data[cmte_key]["pe"] for cmte_key in keys(cmte_data)]

    push!(uq_dict["mean_fcomp_stdevs"],avg_fcomp_stdev)
    push!(uq_dict["energy_stdevs"],energy_stdev)

    if ace_cmte.extra_data 
        push!(uq_dict["each_fcomp_mean"],fcomp_means)
        push!(uq_dict["each_fcomp_stdevs"],fcomp_stdevs)
        push!(uq_dict["mean_pe"],energy_mean)
        push!(uq_dict["cmte_pe_preds"],raw_energies)
    end

    if (avg_fcomp_stdev > ace_cmte.force_thresh) || (energy_stdev > ace_cmte.energy_thresh)
        push!(uncertain_configs,deepcopy(main_cfg_dict))
    end
end

function extended_extract_ss_obs(lmp::LMP)
    atomids    = extract_atom(lmp, "id")
    raw_types  = extract_atom(lmp,"type")
    raw_pos    = extract_atom(lmp, "x")
    raw_forces = extract_atom(lmp, "f")
    raw_vels   = extract_atom(lmp, "v")

    pe  = extract_compute(lmp,"pe",LAMMPS.API.LMP_STYLE_GLOBAL,LAMMPS.API.LMP_TYPE_SCALAR)
    temp  = extract_compute(lmp,"temp",LAMMPS.API.LMP_STYLE_GLOBAL,LAMMPS.API.LMP_TYPE_SCALAR)
    xlo = LAMMPS.extract_variable(lmp, "xlo")
    xhi = LAMMPS.extract_variable(lmp, "xhi")
    ylo = LAMMPS.extract_variable(lmp, "ylo")
    yhi = LAMMPS.extract_variable(lmp, "yhi")
    zlo = LAMMPS.extract_variable(lmp, "zlo")
    zhi = LAMMPS.extract_variable(lmp, "zhi")
    xy = LAMMPS.extract_variable(lmp, "xy")
    xz = LAMMPS.extract_variable(lmp, "xz")
    yz = LAMMPS.extract_variable(lmp, "yz")

    # Enforce orthogonality assumption 
    @assert isapprox(xy,0.0;atol=1e-24)
    @assert isapprox(xz,0.0;atol=1e-24)
    @assert isapprox(yz,0.0;atol=1e-24)


    masses = extract_atom(lmp, "mass")[2:end] # not sure the point of the first entry, maybe a 1-indexing translation thing?
   
    sort_idxs = sortperm(atomids)
    config_types  = raw_types[sort_idxs]
    config_pos    = raw_pos[:, sort_idxs] 
    config_forces = raw_forces[:, sort_idxs]
    config_vels   = raw_vels[:, sort_idxs]

    config_dict = Dict( "types"      => config_types,
                        "positions"  => config_pos,
                        "forces"     => config_forces,
                        "velocities" => config_vels,
                        "box_bounds" => [ [xlo,xhi], [ylo,yhi], [zlo,zhi]],
                        "masses"     => masses,
                        "pe"         => pe,
                        "temp"       => temp,
                        )

    
    config_dict
end

function extract_single_step_observables(lmp::LMP)
    atomids    = extract_atom(lmp, "id")
    raw_types  = extract_atom(lmp,"type")
    raw_pos    = extract_atom(lmp, "x")
    raw_forces = extract_atom(lmp, "f")
    raw_vels   = extract_atom(lmp, "v")

    pe = extract_compute(lmp,"pe",LAMMPS.API.LMP_STYLE_GLOBAL,LAMMPS.API.LMP_TYPE_SCALAR)
    
   
    sort_idxs = sortperm(atomids)
    config_types  = raw_types[sort_idxs]
    config_pos    = raw_pos[:, sort_idxs] 
    config_forces = raw_forces[:, sort_idxs]
    config_vels   = raw_vels[:, sort_idxs]

    config_dict = Dict( "types"      => config_types,
                        "positions"  => config_pos,
                        "forces"     => config_forces,
                        "velocities" => config_vels,
                        "pe"         => pe)
    
    config_dict
end


function ace_singlepoint(bbounds,atom_types,atom_pos,masses, tstep, lmp)
    num_types = size(masses)[1]
    num_atoms = size(atom_pos)[2]

    bbound_str = ""
    cart_name = ["x", "y", "z"]
    for i in 1:3
        bbound_str = bbound_str * " $(cart_name[i]) final" * (@sprintf " %.27f" bbounds[i][1]) * (@sprintf " %.27f" bbounds[i][2])
    end
    bbound_str = bbound_str * " boundary p p p"

    command(lmp, "change_box all" * bbound_str)
    
    for (type_idx,mass) in enumerate(masses)
        command(lmp, "mass $(type_idx) $(mass)")
    end

    # atom_types and atom_pos need to be in the same order, taken care of in extract function
    for j in 1:num_atoms
        xyz_str = ""
        for i in 1:3
            xyz_str = xyz_str * @sprintf " %.27f" atom_pos[i,j]
        end
        command(lmp, "create_atoms $(atom_types[j]) single" * xyz_str)
    end
    
    command(lmp, "reset_timestep $(tstep)")
    command(lmp, "run 0")
    config_dict = extract_single_step_observables(lmp)

    command(lmp, "delete_atoms group all") 
    return config_dict
end



function update_thermo_record!(cfg_dict,thermo_dict)
    push!(thermo_dict["all_temps"],cfg_dict["temp"])
    push!(thermo_dict["all_pes"],cfg_dict["pe"])
end

function md_activelearn(md_initialization::Function,init_param_dict,uq_metric,num_steps=1000)
    lmp, startstop = md_initialization(;init_param_dict...)

    uq_dict = initialize_uq_dict(uq_metric)
    thermo_dict = Dict("all_temps" => [],
                       "all_pes"   => [])
    uncertain_configs = []

    command(lmp, "thermo       1")
    command(lmp, "thermo_style custom step temp pe ke etotal press")

    command(lmp, "reset_timestep 0")

    command(lmp,"compute pe all pe")
    command(lmp,"compute temp all temp")
     
    # implementing get_thermo would obviate the need for this
    # Though maybe I should check extract_global() ???
    command(lmp,"variable xlo equal xlo")
    command(lmp,"variable xhi equal xhi" )
    command(lmp,"variable ylo equal ylo")
    command(lmp,"variable yhi equal yhi" )
    command(lmp,"variable zlo equal zlo")
    command(lmp,"variable zhi equal zhi" )
    command(lmp,"variable xy equal xy")
    command(lmp,"variable xz equal xz")
    command(lmp,"variable yz equal yz")


    command(lmp, "run 0")
    main_cfg_dict = extended_extract_ss_obs(lmp)
    update_thermo_record!(main_cfg_dict,thermo_dict)
    update_uq!(uq_metric,main_cfg_dict,uq_dict,uncertain_configs)

    for tstep in 1:num_steps
        if tstep % 1000 == 0
            println("Timestep: $(tstep)")
        end
        
        try
            run_command = "run 1"
            if !isnothing(startstop)
                run_command = run_command *""" start $(startstop[:start]) stop $(startstop[:stop])"""
            end
            command(lmp, run_command)
        catch e
            break 
        end

        main_cfg_dict = extended_extract_ss_obs(lmp)
        update_thermo_record!(main_cfg_dict,thermo_dict)
        update_uq!(uq_metric,main_cfg_dict,uq_dict,uncertain_configs) 
    end

    # Doing this apparently fucks with the Garbage collection later on 
    #LAMMPS.API.lammps_close(lmp)

    uq_dict,thermo_dict,uncertain_configs
end