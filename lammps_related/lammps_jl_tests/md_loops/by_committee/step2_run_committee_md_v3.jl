using LAMMPS
using CSV, DataFrames
using Plots
using Printf
using Statistics 

function extended_extract_ss_obs(lmp::LMP)
    atomids    = extract_atom(lmp, "id")
    raw_types  = extract_atom(lmp,"type")
    raw_pos    = extract_atom(lmp, "x")
    raw_forces = extract_atom(lmp, "f")
    raw_vels   = extract_atom(lmp, "v")

    pe  = extract_compute(lmp,"pe",LAMMPS.API.LMP_STYLE_GLOBAL,LAMMPS.API.LMP_TYPE_SCALAR)
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
                        "pe"         => pe)
    
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

function initialize_committee_member(ace_fname, num_types)
    lmp = LMP(["-screen", "none"])
    command(lmp, "log none")

    command(lmp, "units          metal")
    command(lmp, "boundary       p p p")
    command(lmp, "atom_style     atomic")
    command(lmp, "neigh_modify   delay 0 every 1 check no") # probably not necessary for a single calc

    command(lmp, "region main block 0.0 1.0 0.0 1.0 0.0 1.0")
    command(lmp, "create_box $(num_types) main")

    command(lmp, "pair_style  pace")
    command(lmp, """pair_coeff     * * $(ace_fname) Hf O""")

    command(lmp, "thermo       1")
    command(lmp, "thermo_style custom step temp pe ke etotal press")
    command(lmp,"compute pe all pe")
    return lmp
end

function compute_avg_force_stdev(cmte_dict)
    f_stdevs = []
    for f_idx in eachindex(cmte_dict[1]["forces"])
        f_comp_std = Statistics.std([cmte_dict[cmte_key]["forces"][f_idx] for cmte_key in keys(cmte_dict)]) # HARDCODED
        push!(f_stdevs, f_comp_std)
    end
    avg_fcomp_stdev = mean(f_stdevs)
    avg_fcomp_stdev
end 

function compute_energy_mean_and_stdev(cmte_dict)
    e_stdev = Statistics.std([ cmte_dict[cmte_key]["pe"] for cmte_key in keys(cmte_dict)])
    e_mean = Statistics.mean([ cmte_dict[cmte_key]["pe"] for cmte_key in keys(cmte_dict)])
    e_mean, e_stdev
end

function ace_committee_expts(driver_yace_fname, cmte_lmps; num_steps=50, vel_seed =12280329, temp=100)

    mean_fcomp_stdevs, energy_stdevs,energy_means = LMP(["-screen","none"]) do lmp  
        data = []
        mean_fcomp_stdevs = []
        energy_stdevs = []
        energy_means = []
        command(lmp, "log none")

        command(lmp, "units          metal")
        command(lmp, "boundary       p p p")
        command(lmp, "atom_style     atomic")
        command(lmp, "neigh_modify   delay 0 every 1 check no") # neighborlist rebuilt every step
        
        command(lmp, "read_data     tetrag_hfo2_sample_DATA")
        
        command(lmp, "pair_style     pace")
        command(lmp, "pair_coeff     * * $(driver_yace_fname) Hf O")
        
        command(lmp, "variable T        equal  $(temp)")
        command(lmp, "variable Tdamp    equal  0.1")
        command(lmp, "variable Tseed    equal  $(vel_seed)")
        command(lmp, "variable dumpf    equal  1")
        
        command(lmp, "velocity     all create \$T \${Tseed} mom yes rot yes dist gaussian")
        
        command(lmp, "thermo       1")
        command(lmp, "thermo_style custom step temp pe ke etotal press")
      
        #command(lmp, "dump           run_forces all custom \${dumpf} dump_CHECK.custom id type x y z fx fy fz vx vy vz")
        #command(lmp, """dump_modify    run_forces sort id format line "%4d %1d %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f" """)
        
        command(lmp,"compute pe all pe")
         
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


        command(lmp, "fix          nvt all nvt temp \$T \$T \${Tdamp}")
        command(lmp, "run 0")

        cmte_data = Dict()
        main_cfg_dct = extended_extract_ss_obs(lmp) 
        cmte_data[1] = main_cfg_dct
        for i in 1:length(cmte_lmps)
            other_cfg_dct = ace_singlepoint(main_cfg_dct["box_bounds"],main_cfg_dct["types"],main_cfg_dct["positions"],main_cfg_dct["masses"],0,cmte_lmps[i])
            cmte_data[i+1] = other_cfg_dct
        end
        avg_fcomp_stdev = compute_avg_force_stdev(cmte_data)
        push!(mean_fcomp_stdevs,avg_fcomp_stdev)
        energy_mean, energy_stdev = compute_energy_mean_and_stdev(cmte_data)
        push!(energy_stdevs,energy_stdev)
        push!(energy_means,energy_mean)

        for tstep in 1:num_steps
            command(lmp, "run 1")
            main_cfg_dct = extended_extract_ss_obs(lmp) 
            cmte_data[1] = deepcopy(main_cfg_dct)
            for i in 1:length(cmte_lmps)
                other_cfg_dct =ace_singlepoint(main_cfg_dct["box_bounds"],main_cfg_dct["types"],main_cfg_dct["positions"],main_cfg_dct["masses"],tstep,cmte_lmps[i])
                cmte_data[i+1] = deepcopy(other_cfg_dct)
            end
            avg_fcomp_stdev = compute_avg_force_stdev(cmte_data)
            push!(mean_fcomp_stdevs,avg_fcomp_stdev)
            energy_mean, energy_stdev = compute_energy_mean_and_stdev(cmte_data)
            push!(energy_stdevs,energy_stdev)
            push!(energy_means,energy_mean)

        end
        return mean_fcomp_stdevs, energy_stdevs, energy_means
    end 
    mean_fcomp_stdevs, energy_stdevs,energy_means
end

label = "HfO2_N2_P6_r4"
yace_files = ["./fits/$(label)_fit$(fit_idx).yace" for fit_idx in 1:5]

yace_lmps = [initialize_committee_member(yace_files[i],2) for i in 2:length(yace_files)] # hard-coding number of types currently...

#@time ace_committee_expts(yace_files[1],yace_lmps; num_steps=5, vel_seed =12280329, temp=100);
# on second run, 1.618354 seconds (43.28 k allocations: 5.352 MiB, 0.55% gc time, 0.19% compilation time)
#@time ace_committee_expts(yace_files[1],yace_lmps; num_steps=50, vel_seed =12280329, temp=100);
# 1.974278 seconds (367.53 k allocations: 45.478 MiB, 0.42% gc time)
#@time ace_committee_expts(yace_files[1],yace_lmps; num_steps=500, vel_seed =12280329, temp=100);
# 5.471430 seconds (3.61 M allocations: 446.726 MiB, 0.66% gc time)

f_stdev_200, e_stdev_200, e_mean_200 = ace_committee_expts(yace_files[1],yace_lmps; num_steps=5000, vel_seed =12280329, temp=200);
f_stdev_500, e_stdev_500, e_mean_500 = ace_committee_expts(yace_files[1],yace_lmps; num_steps=5000, vel_seed =12280329, temp=500);
f_stdev_1000, e_stdev_1000, e_mean_1000 = ace_committee_expts(yace_files[1],yace_lmps; num_steps=5000, vel_seed =12280329, temp=1000);
f_stdev_2000, e_stdev_2000, e_mean_2000 = ace_committee_expts(yace_files[1],yace_lmps; num_steps=5000, vel_seed =12280329, temp=2000);

plot(1:5001, [f_stdev_200,f_stdev_500,f_stdev_1000,f_stdev_2000])
plot(1:5001, [e_stdev_200,e_stdev_500,e_stdev_1000,e_stdev_2000])
plot(1:5001, [e_mean_200,e_mean_500,e_mean_1000,e_mean_2000])

for c_lmps in yace_lmps
    LAMMPS.API.lammps_close(c_lmps)
end

