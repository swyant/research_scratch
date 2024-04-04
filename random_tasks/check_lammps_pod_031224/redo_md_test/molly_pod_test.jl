using AtomsIO 
using Molly 
using Atomistic 
using InteratomicPotentials
using GLMakie

spot_check_sys = load_system("./md_test/spot_check_molly80.xyz")

hfo2_sys_raw = load_system("./md_test/sample_monoclinic_HfO2.xyz")
hfo2_sys = staticAtoms(hfo2_sys_raw)
molly_hfo2_sys = Molly.System(hfo2_sys, u"eV", u"eV/Å")

lmp_pod1 = LAMMPS_POD("./sample_6body_hfo2_param.pod", [:Hf,:O])
hfo2_lbp1 = LBasisPotential(lmp_pod1,"./sample_6body_2elem_coeffs.pod")

lmp_pod2 = LAMMPS_POD("./sample_6body_hfo2_param.pod", [:Hf,:O])
hfo2_lbp2 = LBasisPotential(lmp_pod2,"./sample_6body_2elem_coeffs.pod")

pod_inter = InteratomicPotentialInter(hfo2_lbp1, u"eV", u"Å")
general_inters = (pod_inter,)

m_sys = Molly.System(molly_hfo2_sys;
                    general_inters=general_inters, 
                    loggers=(#force=ForceLogger(typeof(1.0u"eV/Å"), 1),
                             #energy=PotentialEnergyLogger(typeof(1.0u"eV"),1),
                             coords=CoordinateLogger(1),
                             vels=VelocityLogger(1))
                    )
         
#simulator = VelocityVerlet(
#    dt=0.001u"ps",
#    remove_CM_motion=0,
#)

#simulator = VelocityVerlet(
#    dt=0.001u"ps",
#    coupling=AndersenThermostat(300u"K", 0.1u"ps"),
#    remove_CM_motion=0,
#)

#random_velocities!(m_sys,200.0u"K")
random_velocities!(m_sys,500.0u"K")

#simulator = Langevin(
#    dt=0.001u"ps",
#    temperature=200u"K",
#    friction=1.0u"ps^-1",
#    remove_CM_motion=0,
#)

simulator = NoseHoover(
    dt=0.001u"ps",
    #temperature=200u"K",
    temperature=500u"K",
    #remove_CM_motion=0
    remove_CM_motion=1
)

color_map = Dict(:Hf => :grey, :O => :red)
size_map = Dict(:Hf => 0.5, :O => 0.25)
colors = [color_map[sym] for sym in atomic_symbol(m_sys)]
sizes = [size_map[sym] for sym in atomic_symbol(m_sys)]

#simulate!(m_sys, simulator,1000)
#visualize(m_sys.loggers.coords, m_sys.boundary, "md_test/test.mp4"; color=colors, markersize=sizes)
