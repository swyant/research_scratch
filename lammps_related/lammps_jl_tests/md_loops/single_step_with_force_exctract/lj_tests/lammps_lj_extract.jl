using LAMMPS
using JLD

single = true
num_steps = 50

LAMMPS.locate()
lmp = LMP(["-screen","none"])

command(lmp, "units          metal")
command(lmp, "boundary       p p p")
command(lmp, "atom_style     atomic")
command(lmp, "neigh_modify   delay 0 every 1 check no")

command(lmp, "read_data     fcc_lj_Ar_smaller_DATA")

command(lmp, "pair_style     lj/cut 8.51")
command(lmp, "pair_coeff     * * 0.01032 3.405")

command(lmp, "variable T        equal  100")
command(lmp, "variable Tdamp    equal  0.1")
command(lmp, "variable Tseed    equal  12280329")
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

command(lmp,"run 1")
raw_forces = extract_atom(lmp, "f")
raw_pos    = extract_atom(lmp, "x")
raw_types  = extract_atom(lmp,"type")
atomids = extract_atom(lmp, "id")

forces = raw_forces[:,sortperm(atomids)]
pos    = raw_pos[:,sortperm(atomids)]
types  = raw_types[sortperm(atomids)]

JLD.save("step1_ordered_extracts.jld", "forces", forces, "pos", pos, "types", types)

if single
    for i in 1:num_steps
        command(lmp, "run 1")
    end
else
    command(lmp, "run $(num_steps)")
end


#id2idx(id::Integer, ids::Vector{<:Integer}) = findfirst(x -> x==id,ids)

