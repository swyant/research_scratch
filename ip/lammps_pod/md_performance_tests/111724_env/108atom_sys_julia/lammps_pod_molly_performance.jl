using AtomsIO 
using Molly 
using Atomistic 
using InteratomicPotentials
#using GLMakie

temperature = 500u"K"

hfo2_sys = load_system("../../files/HfO2_108atom.xyz")
molly_hfo2_sys = Molly.System(hfo2_sys; energy_units=u"eV", force_units=u"eV/Å")

lmp_pod = LAMMPS_POD("../../files/sample_6body_hfo2_param.pod", [:Hf,:O])
hfo2_lbp = LBasisPotential(lmp_pod,"../../files/sample_6body_2elem_coeffs.pod")

pod_inter = InteratomicPotentialInter(hfo2_lbp, u"eV", u"Å")
general_inters = (pod_inter,)

m_sys = Molly.System(molly_hfo2_sys;
                    general_inters=general_inters, 
                    )
         
simulator = NoseHoover(
    dt=0.001u"ps",
    temperature=temperature,
    remove_CM_motion=1
)

simulate!(m_sys, simulator,1)

for _ in 1:5
  run_sys = deepcopy(m_sys)
  random_velocities!(run_sys,temperature)
  
  @time simulate!(run_sys, simulator,3000)
  @show run_sys.coords[1:20]
end