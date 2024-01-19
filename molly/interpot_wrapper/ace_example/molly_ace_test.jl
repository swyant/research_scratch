using CSV, DataFrames
cd("/Users/swyant/cesmix/exploratory/new_public/molly/interpot_wrapper/ace_example")

include("../molly_interpot_utils.jl")

ace = ACE(  species            = [:Hf, :O],
            body_order         = 4, 
            polynomial_degree  = 10,
            wL                 = 1.5, 
            csp                = 1.0,
            r0                 = 2.15,
            rcutoff            = 5.0)

coeffs = vec(Matrix(CSV.read("./N3_rcut5_maxdeg10_1e-3lambdaQR_fit_coeffs.csv",DataFrame,header=false)))

lb = LBasisPotential(coeffs,[0.0],ace)

inter_ace = InteratomicPotentialInter(lb,InteratomicPotentials.energy_and_force)
general_inters = (inter_ace,)

sys = load_system("test_tetrag_hfo2.xyz")

compute_local_descriptors(sys,ace)
compute_force_descriptors(sys,ace)
InteratomicPotentials.force(sys, lb)

f_check = Molly.forces(inter_ace, sys)

mp = molly_params(sys)
m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(typeof(NoUnits), 1),),
            force_units=NoUnits,
            energy_units=NoUnits,
            #loggers=(force=ForceLogger(Float32, 1),),
            #force_units=NoUnits,
            )

simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(m_sys,simulator,100)

