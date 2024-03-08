using LAMMPS
lmp = LMP(["-screen", "none"])


command(lmp, "units          metal")
command(lmp, "boundary       p p p") 
command(lmp, "atom_style     atomic")

#command(lmp, "read_data      ../test_monoclinic_hfo2_unit_DATA")

#This is causing the error (instead of reading the DATA fileit)
# Also changed the name from main to dummy
command(lmp, "region dummy prism 0.0 2.0 0.0 2.0 0.0 2.0 0.0 0.0 0.0")
command(lmp, "create_box 2 dummy")

#This made it work
#command(lmp, "region dummy prism 0.0 5.1364755065543157 0.0 5.1934133375444702 0.0 5.2461598045212767 0.0 -0.88344296827430302 0.0")
#command(lmp, "create_box 2 dummy")

#command(lmp, "region dummy prism 0.0 5.1364755065543157 0.0 5.1934133375444702 0.0 5.2461598045212767 0.0 0.0 0.0")
#command(lmp, "create_box 2 dummy")


command(lmp, "mass 1 178.49")
command(lmp, "mass 2 15.999")

command(lmp, "create_atoms 1 single 0.5 0.5 0.5")


# This made it work
#command(lmp, "create_atoms 1 single 0.79296805999999997     2.3789195900000002      3.7137930099999998")
#command(lmp, "create_atoms 1 single 3.0183409399999999      4.9756267799999998      4.1554425899999998")
#command(lmp, "create_atoms 1 single 3.4600634499999998              2.81449477      1.5323642200000001")
#command(lmp, "create_atoms 1 single 1.2346900599999999              0.21778759              1.09071465")
#command(lmp, "create_atoms 2 single   1.4366105600000001      3.8562449399999998      5.1337742300000002")
#command(lmp, "create_atoms 2 single           2.37469844              1.25953776      2.7354618799999999")
#command(lmp, "create_atoms 2 single   2.8164204399999999      1.3371694199999999                0.112383")
#command(lmp, "create_atoms 2 single           1.87833256              3.93387661      2.5106953600000002")
#command(lmp, "create_atoms 2 single 0.045238479999999998              1.71776788      1.8164965900000001")
#command(lmp, "create_atoms 2 single   4.6495150299999999      4.3144750600000004     0.80658229000000004")
#command(lmp, "create_atoms 2 single   4.2077925199999999      3.4756464899999999      3.4296606500000002")
#command(lmp, "create_atoms 2 single -0.39648402999999999     0.87893931000000003      4.4395749499999999")

#command(lmp, "create_atoms 2 single 0.75 0.75 0.75")

command(lmp, "neighbor       0.5 bin") 
command(lmp, "neigh_modify  delay 0 every 1 check no")    


#command(lmp, "pair_style     pod")
#command(lmp, """pair_coeff     * * ../HfO2_FPOD_020224_param.pod ../HfO2_FPOD_020224_v2_coefficients.pod "" "" Hf O""")
command(lmp, "pair_style    zero 10.0")
command(lmp, "pair_coeff    * * ")


command(lmp, """compute dd all podd/atom ../HfO2_FPOD_020224_param.pod "" "" Hf O""")

command(lmp, "dump 		     mydump_dd all custom 1 dump_dd id type c_dd[*]")
command(lmp, "dump_modify  mydump_dd sort id format float %20.15g")

command(lmp, "run 0")
