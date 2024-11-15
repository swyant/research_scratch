using AtomsIO, Molly  # JuliaMolSim packages
using Atomistic, InteratomicPotentials # CESMIX packages

hfo2_sys = load_system("./HfO2_108atom.xyz")
molly_hfo2_sys = Molly.System(hfo2_sys; 
                              energy_units=u"eV", 
                              force_units=u"eV/Å")

# construct (LAMMPS-implemented) POD potential from IP.jl
lmp_pod = LAMMPS_POD("./sample_6body_hfo2_param.pod", [:Hf,:O])
hfo2_lbp = LBasisPotential(lmp_pod,"./sample_6body_2elem_coeffs.pod")

# Wrapper compatible with AtomsCalculators.jl
pod_inter = InteratomicPotentialInter(hfo2_lbp, u"eV", u"Å")
general_inters = (pod_inter,)

# Example of directly using AtomsCalculators.jl interface
forces = AtomsCalculators.forces(molly_hfo2_sys,pod_inter)
pe     = AtomsCalculators.potential_energy(molly_hfo2_sys,pod_inter)

# Set up Molly and run short MD
m_sys = Molly.System(molly_hfo2_sys;
                    general_inters=general_inters, 
                    )

random_velocities!(run_sys,temperature)
simulator = NoseHoover(
    dt=0.001u"ps",
    temperature=temperature,
    remove_CM_motion=1
)

simulate!(m_sys, simulator,10_000)
















