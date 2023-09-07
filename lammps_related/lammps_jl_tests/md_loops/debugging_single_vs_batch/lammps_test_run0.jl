using LAMMPS

LAMMPS.locate()
lmp = LMP(["-screen","none"])

command(lmp, "units          metal")
command(lmp, "boundary       p p p")
command(lmp, "atom_style     atomic")

command(lmp, "read_data     discrepant_hfo2_config_DATA")

command(lmp, "pair_style     pace")
command(lmp, "pair_coeff     * * Siviraman_GAP_HfO2_N3_rcut5_maxdeg10_regQR.yace Hf O")

command(lmp, "dump           final_dump all custom 1 dump_run0.custom id type x y z fx fy fz")
command(lmp, """dump_modify  final_dump sort id format line "%4d %1d %21.16f %21.16f %21.16f %21.16f %21.16f %21.16f" """)

command(lmp, "run 0")
