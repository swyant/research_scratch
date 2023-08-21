using LAMMPS 

# Just do this once, then restart
#LAMMPS.set_library!("/Users/swyant/local_software/lammps/lammps_081723/shared_build/liblammps.0.dylib")

LAMMPS.locate()
#lmp = LMP(["-screen","none"])
#
#command(lmp, "units  metal")
#command(lmp, "boundary  p p p")
#command(lmp, "atom_style  atomic")
#
#command(lmp, "read_data  tetrag_hfo2_sample_DATA")
#
#command(lmp, "pair_style  quip")
#command(lmp, """pair_coeff  * * gap.xml "Potential xml_label=GAP_2020_2_11_0_18_44_47_601" 72 8""")
#
#command(lmp, "variable T  equal  2000")
#command(lmp, "variable Tdamp  equal  0.1")
#command(lmp, "variable Tseed  equal  12280329")
#command(lmp, "variable dumpf  equal  50")
#
#command(lmp, "velocity  all create \$T \${Tseed} mom yes rot yes dist gaussian")
#
#command(lmp, "thermo  10")
#command(lmp, "thermo_style custom step temp pe ke etotal press")
#
#command(lmp, "dump  run_dump all xyz \${dumpf} dump_run.xyz")
#command(lmp, """dump_modify  run_dump sort id format line "%s %21.16f %21.16f %21.16f" """)
#
#command(lmp, "dump  run_forces all custom \${dumpf} dump_run_forces.custom id type x y z fx fy fz")
#command(lmp, """dump_modify    run_forces sort id format line "%4d %1d %21.16f %21.16f %21.16f %21.16f %21.16f %21.16f" """)
#
#
#command(lmp, "fix  nvt all nvt temp \$T \$T ${Tdamp}")
#command(lmp, "run 500")