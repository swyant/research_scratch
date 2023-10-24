using LAMMPS
using CSV, DataFrames
#using JLD

include("../../../utilities/parse_dump/parse_dump.jl")
include("../lammps_extract_utils.jl")

function lj_expts(; num_steps=50, vel_seed =12280329, single=true)
    configs = LMP(["-screen","none"]) do lmp  
        command(lmp, "units          metal")
        command(lmp, "boundary       p p p")
        command(lmp, "atom_style     atomic")
        command(lmp, "neigh_modify   delay 0 every 1 check no") # neighborlist rebuilt every step
        
        command(lmp, "read_data     fcc_lj_Ar_smaller_DATA")
        
        command(lmp, "pair_style     lj/cut 8.51")
        command(lmp, "pair_coeff     * * 0.01032 3.405")
        
        command(lmp, "variable Tstart        equal  10")
        command(lmp, "variable Tend          equal 200")
        command(lmp, "variable Tdamp    equal  0.1")
        command(lmp, "variable Tseed    equal  $(vel_seed)")
        command(lmp, "variable dumpf    equal  1")
        
        command(lmp, "velocity     all create \${Tstart} \${Tseed} mom yes rot yes dist gaussian")
        
        command(lmp, "thermo       1")
        command(lmp, "thermo_style custom step temp pe ke etotal press")
        
        if single
            command(lmp, "dump           run_forces all custom \${dumpf} dump_single_nvt_ramp.custom id type x y z fx fy fz vx vy vz")
        else
            command(lmp, "dump           run_forces all custom \${dumpf} dump_batch_nvt_ramp.custom id type x y z fx fy fz vx vy vz")
        end
        command(lmp, """dump_modify    run_forces sort id format line "%4d %1d %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f" """)
        
        if single
            command(lmp,"compute pe all pe")
        else
            command(lmp,"""fix eprint all print 1 "\$(step) \$(pe)" file pe.dat""")
        end
        command(lmp, "fix          nvt all nvt temp \${Tstart} \${Tend} \${Tdamp}")
         
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
            
            configs = parse_pe_file(;config_list=configs)
            
            #rm("./dump_batch.custom") # could also just be a file in /tmp
        end
        return configs
    end 
    configs
end


single_step_configs = lj_expts(; num_steps=5000, single=true)
ss_all_forces = reshape(reduce(hcat,[config["forces"] for config in single_step_configs]), size(single_step_configs[1]["forces"])..., :) 
ss_all_pe = [config["pe"] for config in single_step_configs]

batch_step_configs = lj_expts(; num_steps=5000, single=false)
batch_all_forces = reshape(reduce(hcat,[config["forces"] for config in batch_step_configs]), size(batch_step_configs[1]["forces"])..., :)
batch_all_pe = [config["pe"] for config in batch_step_configs]

@show max_force_discrep = maximum(abs.(batch_all_forces - ss_all_forces)) # 4.985600920996795e-28
@show max_pe_discrep = maximum(abs.(batch_all_pe - ss_all_pe)) # 0.0