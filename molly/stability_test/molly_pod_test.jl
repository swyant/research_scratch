using AtomsIO 
using Molly 
using Atomistic 
using InteratomicPotentials
using CSV, DataFrames
using GLMakie

#hfo2_sys_raw = load_system("./sample_monoclinic_HfO2.xyz")
hfo2_sys_raw = load_system("./sample_monoclinic_HfO2_wdefect.xyz")

hfo2_sys = staticAtoms(hfo2_sys_raw)
molly_hfo2_sys = Molly.System(hfo2_sys, u"eV", u"eV/Å")

#lmp_pod = LAMMPS_POD("./HfO2_FPOD_020224_param.pod", [:Hf,:O])
#hfo2_lbp = LBasisPotential(lmp_pod,"./HfO2FPOD2_coefficients.pod")

ace_hfo2 = ACE(species            = [:Hf, :O],
               body_order         = 4, 
               polynomial_degree  = 10,
               wL                 = 1.5, 
               csp                = 1.0,
               r0                 = 2.15,
               rcutoff            = 5.0)

coeffs_hfo2 = vec(Matrix(CSV.read("./N3_rcut5_maxdeg10_1e-3lambdaQR_fit_coeffs.csv",DataFrame,header=false)))

hfo2_lbp = LBasisPotential(coeffs_hfo2,[0.0],ace_hfo2)


pod_inter = InteratomicPotentialInter(hfo2_lbp, u"eV", u"Å")
general_inters = (pod_inter,)

m_sys = Molly.System(molly_hfo2_sys;
                    general_inters=general_inters, 
                    loggers=(#force=ForceLogger(typeof(1.0u"eV/Å"), 1),
                             #energy=PotentialEnergyLogger(typeof(1.0u"eV"),1),
                             coords=CoordinateLogger(1),
                             vels=VelocityLogger(1))
                    )
         
random_velocities!(m_sys,3000.0u"K")

simulator = NoseHoover(
    dt=0.001u"ps",
    temperature=500u"K",
    remove_CM_motion=1
)

color_map = Dict(:Hf => :grey, :O => :red)
size_map = Dict(:Hf => 0.5, :O => 0.25)
colors = [color_map[sym] for sym in atomic_symbol(m_sys)]
sizes = [size_map[sym] for sym in atomic_symbol(m_sys)]


#simulate!(m_sys, simulator,1000)
#visualize(m_sys.loggers.coords, m_sys.boundary, "./test.mp4"; color=colors, markersize=sizes)
