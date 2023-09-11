using LAMMPS

LAMMPS.locate()
lmp = LMP(["-screen","none"])

command(lmp, "units          metal")
command(lmp, "boundary       p p p")
command(lmp, "atom_style     atomic")
command(lmp, "neigh_modify   delay 0 every 1 check no")

command(lmp, "read_data     fcc_lj_Ar_box-zeroed_shifted_DATA")

command(lmp, "pair_style     lj/cut 8.51")
command(lmp, "pair_coeff     * * 0.01032 3.405")

command(lmp, "variable T        equal  2000")
command(lmp, "variable Tdamp    equal  0.1")
command(lmp, "variable Tseed    equal  12280329")
command(lmp, "variable dumpf    equal  1")

command(lmp, "velocity     all create \$T \${Tseed} mom yes rot yes dist gaussian")

command(lmp, "thermo       1")
command(lmp, "thermo_style custom step temp pe ke etotal press")

command(lmp, "dump           run_forces all custom \${dumpf} dump_single.custom id type x y z fx fy fz vx vy vz")
command(lmp, """dump_modify    run_forces sort id format line "%4d %1d %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f" """)


command(lmp, "fix          nvt all nvt temp \$T \$T \${Tdamp}")

for i in 1:50
    command(lmp, "run 1")
end
