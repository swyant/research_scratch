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
   command(lmp, "neigh_modify   delay 0 every 1 check no")

   command(lmp, "read_data    tetrag_hfo2_sample_DATA")
   command(lmp, "pair_style   pace")
   command(lmp, "pair_coeff   * * $(driver_yace_fname)$(atomtype_str)")
   
   command(lmp, "velocity     all create $(start_temp) $(vel_seed)"*
                              " mom yes rot yes dist gaussian")
   command(lmp, "fix          equib_npt all npt temp $(start_temp)"*
                              " $(start_temp) 0.1 iso 0.0 0.0 1.0")
   command(lmp, "run          $(equib_steps)")   
   command(lmp, "unfix        equib_npt")
   command(lmp, "reset_timestep 0")
 
   command(lmp, "fix          npt_ramp all npt temp $(start_temp)"*
                              " $(end_temp) 0.1 iso 0.0 0.0 1.0")
   start_stop_dict = Dict(:start => "0",
                          :stop  => "$(ramp_steps)")
   lmp, start_stop_dict
end

yace_files = ["$(fit_dir)/$(label)_fit$(fit_idx).yace" for fit_idx in 1:5]
ace_cmte = PACE_Ensemble(yace_files, ["Hf","O"]; energy_thresh=0.1)

npt_melt_params = Dict(:driver_yace_fname => yace_files[1],
                      :start_temp         => 2800,
                      :end_temp           => 4000,
                      :vel_seed           => 12280329,
                      :ramp_steps         => 200000,
                      :atom_type_list    => ["Hf","O"],
                      )

uq, thermo, configs = md_activelearn(simple_ace_iso_npt_melt,
                                     npt_melt_params,
                                     ace_cmte,
                                     200000) 
