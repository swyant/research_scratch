using LAMMPS
using CSV, DataFrames
#using JLD

include("../../utilities/parse_dump/parse_dump.jl")
include("./single_point_gap.jl")

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


function ace_and_gap_expts(; num_steps=50, vel_seed =12280329, temp=100)
    ace_configs, gap_configs = LMP(["-screen","none"]) do lmp  
        command(lmp, "log none")

        command(lmp, "units          metal")
        command(lmp, "boundary       p p p")
        command(lmp, "atom_style     atomic")
        command(lmp, "neigh_modify   delay 0 every 1 check no") # neighborlist rebuilt every step
        
        command(lmp, "read_data     tetrag_hfo2_sample_DATA")
        
        command(lmp, "pair_style     pace")
        command(lmp, "pair_coeff     * * Siviraman_GAP_HfO2_N3_rcut5_maxdeg10_regQR.yace Hf O")
        
        command(lmp, "variable T        equal  $(temp)")
        command(lmp, "variable Tdamp    equal  0.1")
        command(lmp, "variable Tseed    equal  $(vel_seed)")
        command(lmp, "variable dumpf    equal  1")
        
        command(lmp, "velocity     all create \$T \${Tseed} mom yes rot yes dist gaussian")
        
        command(lmp, "thermo       1")
        command(lmp, "thermo_style custom step temp pe ke etotal press")
        
        command(lmp, "dump           run_forces all custom \${dumpf} dump_single.custom id type x y z fx fy fz vx vy vz")
        command(lmp, """dump_modify    run_forces sort id format line "%4d %1d %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f" """)
        
        command(lmp,"compute pe all pe")
         
        # implementing get_thermo would obviate the need for this
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
         
        ace_configs = []
        gap_configs = []
        command(lmp, "run 0")
        cfg_dct = extended_extract_ss_obs(lmp)
        gap_cfg_dct = get_gap_eandf(cfg_dct["box_bounds"],cfg_dct["types"],cfg_dct["positions"],cfg_dct["masses"])

        push!(ace_configs,cfg_dct)
        push!(gap_configs,gap_cfg_dct)

        for _ in 1:num_steps
            command(lmp, "run 1")
            cfg_dct = extended_extract_ss_obs(lmp)
            gap_cfg_dct = get_gap_eandf(cfg_dct["box_bounds"],cfg_dct["types"],cfg_dct["positions"],cfg_dct["masses"])

            push!(ace_configs, cfg_dct)
            push!(gap_configs,gap_cfg_dct)
        end
        return ace_configs, gap_configs
    end 
    ace_configs, gap_configs
end


ace_cfgs, gap_cfgs = ace_and_gap_expts(; num_steps=1, temp=100)

# need to recheck how I fit ACE, because I thought I fit to cohesive energy?
ace_pes = [config["pe"] for config in ace_cfgs]
gap_pes = [config["pe"] for config in gap_cfgs] 