using LAMMPS

LAMMPS.locate()
lmp = LMP(["-screen","none"])

command(lmp, "units          metal")
command(lmp, "boundary       p p p")
command(lmp, "atom_style     atomic")
command(lmp, "neigh_modify   delay 0 every 1 check no")

command(lmp, "read_data     tetrag_hfo2_sample_DATA")

command(lmp, "pair_style     pace")
command(lmp, "pair_coeff     * * Siviraman_GAP_HfO2_N3_rcut5_maxdeg10_regQR.yace Hf O")

command(lmp, "variable T        equal  2000")
command(lmp, "variable Tdamp    equal  0.1")
command(lmp, "variable Tseed    equal  12280329")
command(lmp, "variable dumpf    equal  1")

command(lmp, "velocity     all create \$T \${Tseed} mom yes rot yes dist gaussian")

command(lmp, "thermo       1")
command(lmp, "thermo_style custom step temp pe ke etotal press")

command(lmp, "dump           run_forces all custom \${dumpf} dump_batch.custom id type x y z fx fy fz vx vy vz")
command(lmp, """dump_modify    run_forces sort id format line "%4d %1d %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f" """)


command(lmp, "fix          nvt all nvt temp \$T \$T \${Tdamp}")
command(lmp, "run 3")
