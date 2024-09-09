using AtomsIO 
using Molly 
using Atomistic 
using InteratomicPotentials
using GLMakie

temperature = 500u"K"

hfo2_sys_raw = load_system("../../files/sample_monoclinic_HfO2.xyz")
hfo2_sys = staticAtoms(hfo2_sys_raw)
molly_hfo2_sys = Molly.System(hfo2_sys, u"eV", u"eV/Å")

lmp_pod = LAMMPS_POD("../../files/sample_6body_hfo2_param.pod", [:Hf,:O])
hfo2_lbp = LBasisPotential(lmp_pod,"../../files/sample_6body_2elem_coeffs.pod")

pod_inter = InteratomicPotentialInter(hfo2_lbp, u"eV", u"Å")
general_inters = (pod_inter,)

m_sys = Molly.System(molly_hfo2_sys;
                    general_inters=general_inters, 
                    loggers=(coords=CoordinateLogger(1),
                             vels=VelocityLogger(1))
                    )
         
random_velocities!(m_sys,temperature)

simulator = NoseHoover(
    dt=0.001u"ps",
    temperature=temperature,
    remove_CM_motion=1
)

color_map = Dict(:Hf => :grey, :O => :red)
size_map = Dict(:Hf => 0.5, :O => 0.25)
colors = [color_map[sym] for sym in atomic_symbol(m_sys)]
sizes = [size_map[sym] for sym in atomic_symbol(m_sys)]

#simulate!(m_sys, simulator,3000)
#visualize(m_sys.loggers.coords, m_sys.boundary, "check.mp4"; color=colors, markersize=sizes)
