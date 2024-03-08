using LAMMPS
lmp = LMP(["-screen", "none"])


command(lmp, "units          metal")
command(lmp, "boundary       p p p") 
command(lmp, "atom_style     atomic")

command(lmp, "neighbor       0.5 bin") 
command(lmp, "neigh_modify  delay 0 every 1 check no")    


command(lmp, "region main prism 0.0 2.0 0.0 2.0 0.0 2.0 0.0 0.0 0.0")
command(lmp, "create_box 2 main")

command(lmp, "mass 1 178.49")
command(lmp, "mass 2 15.999")

command(lmp, "create_atoms 1 single 0.5 0.5 0.5")

command(lmp, "pair_style    zero 10.0")
command(lmp, "pair_coeff    * * ")


command(lmp, """compute dd all podd/atom ../HfO2_FPOD_020224_param.pod "" "" Hf O""")

command(lmp, "dump 		     mydump_dd all custom 1 dump_dd id type c_dd[*]")
command(lmp, "dump_modify  mydump_dd sort id format float %20.15g")

#command(lmp, "run 0")
