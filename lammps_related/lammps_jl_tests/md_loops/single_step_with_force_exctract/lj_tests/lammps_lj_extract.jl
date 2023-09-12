using LAMMPS
#using JLD

include("../../../utilities/parse_dump/parse_dump.jl")

function extract_single_step_observables(lmp::LMP)
    atomids    = extract_atom(lmp, "id")
    raw_types  = extract_atom(lmp,"type")
    raw_pos    = extract_atom(lmp, "x")
    raw_forces = extract_atom(lmp, "f")
    raw_vels   = extract_atom(lmp, "v")
   
    sort_idxs = sortperm(atomids)
    config_types  = raw_types[sort_idxs]
    config_pos    = raw_pos[:, sort_idxs] 
    config_forces = raw_forces[:, sort_idxs]
    config_vels   = raw_vels[:, sort_idxs]

    config_dict = Dict( "types"      => config_types,
                        "positions"  => config_pos,
                        "forces"     => config_forces,
                        "velocities" => config_vels)
    
    config_dict
end

function lj_expts(; num_steps=50, vel_seed =12280329, single=true)
    configs = LMP(["-screen","none"]) do lmp  
        command(lmp, "units          metal")
        command(lmp, "boundary       p p p")
        command(lmp, "atom_style     atomic")
        command(lmp, "neigh_modify   delay 0 every 1 check no") # neighborlist rebuilt every step
        
        command(lmp, "read_data     fcc_lj_Ar_smaller_DATA")
        
        command(lmp, "pair_style     lj/cut 8.51")
        command(lmp, "pair_coeff     * * 0.01032 3.405")
        
        command(lmp, "variable T        equal  100")
        command(lmp, "variable Tdamp    equal  0.1")
        command(lmp, "variable Tseed    equal  $(vel_seed)")
        command(lmp, "variable dumpf    equal  1")
        
        command(lmp, "velocity     all create \$T \${Tseed} mom yes rot yes dist gaussian")
        
        command(lmp, "thermo       1")
        command(lmp, "thermo_style custom step temp pe ke etotal press")
        
        if single
            command(lmp, "dump           run_forces all custom \${dumpf} dump_single.custom id type x y z fx fy fz vx vy vz")
        else
            command(lmp, "dump           run_forces all custom \${dumpf} dump_batch.custom id type x y z fx fy fz vx vy vz")
        end
        command(lmp, """dump_modify    run_forces sort id format line "%4d %1d %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f" """)
            
        command(lmp, "fix          nvt all nvt temp \$T \$T \${Tdamp}")
         
        configs = []
        if single
            command(lmp, "run 0")
            config_dict = extract_single_step_observables(lmp)
            push!(configs,config_dict)

            for _ in 1:num_steps
                command(lmp, "run 1")
                config_dict = extract_single_step_observables(lmp)
                push!(configs, config_dict)
            end

        else
            command(lmp, "run $(num_steps)")
            configs = fragile_parse_dump("./dump_batch.custom")
            #rm("./dump_batch.custom") # could also just be a file in /tmp
        end
        return configs
    end 
    configs
end


#id2idx(id::Integer, ids::Vector{<:Integer}) = findfirst(x -> x==id,ids)
#cd("/Users/swyant/cesmix/exploratory/new_public/lammps_related/lammps_jl_tests/md_loops/single_step_with_force_exctract/lj_tests")
single_step_configs = lj_expts(; num_steps=1, single=true)
ss_all_forces = reshape(reduce(hcat,[config["forces"] for config in single_step_configs]), size(single_step_configs[1]["forces"])..., :) 

batch_step_configs = lj_expts(; num_steps=1, single=false)
batch_all_forces = reshape(reduce(hcat,[config["forces"] for config in batch_step_configs]), size(batch_step_configs[1]["forces"])..., :)

max_force_discrep = maximum(abs.(batch_all_forces - ss_all_forces))