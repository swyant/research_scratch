# load things 
include("./new_code/dependencies.jl")
include("./new_code/md_active_learning.jl")
include("./new_code/train_utils.jl")

function simple_ace_iso_npt_melt(;driver_yace_fname="sample.yace",
    start_temp=300, 
    end_temp=36000,
    equib_steps = 10000,
    ramp_steps = 200000,
    vel_seed=12280329, 
    atom_type_list=["Hf","O"])

atomtype_str = ""
for elem_str in atom_type_list
atomtype_str = atomtype_str * " " * elem_str
end


lmp = LMP(["-screen","none"])
command(lmp, "log none")

command(lmp, "units          metal")
command(lmp, "boundary       p p p")
command(lmp, "atom_style     atomic")
command(lmp, "neigh_modify   delay 0 every 1 check no") # neighborlist rebuilt every step

command(lmp, "read_data     tetrag_hfo2_sample_DATA")

command(lmp, "pair_style     pace")
command(lmp, "pair_coeff     * * $(driver_yace_fname)$(atomtype_str)")

#command(lmp, "variable Tdamp    equal  0.1")
#command(lmp, "variable dumpf    equal  100")

command(lmp, "velocity     all create $(start_temp) $(vel_seed) mom yes rot yes dist gaussian")
command(lmp, "fix          equib_npt all npt temp  $(start_temp) $(start_temp) 0.1 iso 0.0 0.0 1.0")
command(lmp, "run          $(equib_steps)")

command(lmp, "unfix        equib_npt")
command(lmp, "reset_timestep 0")

#command(lmp, "dump             run_forces all custom 1 dump_unstable_melt.custom id type x y z fx fy fz")
#command(lmp, """dump_modify    run_forces append yes sort id format line "%4d %1d %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f" """)

command(lmp, "fix          npt_ramp all npt temp $(start_temp) $(end_temp) 0.1 iso 0.0 0.0 1.0")

println("successfully ran the initialization commands")

start_stop_dict = Dict(:start => "0",
:stop  => "$(ramp_steps)")

lmp, start_stop_dict
end  


fit_dir = "./initial_fits"
label = "HfO2_N2_P6_r4"

yace_files = ["$(fit_dir)/$(label)_fit$(fit_idx).yace" for fit_idx in 1:5]
ace_cmte = PACE_Ensemble(yace_files, ["Hf","O"]; energy_thresh=0.1,extra_data=true)

npt_melt_params = Dict(:driver_yace_fname => yace_files[1],
                       :start_temp         => 2800,
                       :end_temp           => 4000,
                       :vel_seed           => 12280329,
                       :ramp_steps         => 200000,
                       :atom_type_list    => ["Hf","O"],
                       )

uq_data, thermo_data, uq_configs = md_activelearn(simple_ace_iso_npt_melt,npt_melt_params,ace_cmte,200000) 

#CSV.write("unstable_fcomp_stdevs.csv",DataFrame( Tables.table(reduce(vcat,uq_data["each_fcomp_stdevs"])) ) ,header=false, delim=" ")
#CSV.write("unstable_fcomp_mean.csv",DataFrame( Tables.table(reduce(vcat,uq_data["each_fcomp_mean"])) ) ,header=false, delim=" ")