using CSV, DataFrames
using AtomsIO
cd("/Users/swyant/cesmix/exploratory/new_public/molly/interpot_wrapper/ace_example")

include("../v2_molly_interpot_utils.jl")

ace = ACE(  species            = [:Hf, :O],
            body_order         = 4, 
            polynomial_degree  = 10,
            wL                 = 1.5, 
            csp                = 1.0,
            r0                 = 2.15,
            rcutoff            = 5.0)

coeffs = vec(Matrix(CSV.read("./N3_rcut5_maxdeg10_1e-3lambdaQR_fit_coeffs.csv",DataFrame,header=false)))

lb = LBasisPotential(coeffs,[0.0],ace)

inter_ace = InteratomicPotentialInter(lb,u"eV", u"Å")
general_inters = (inter_ace,)

sys = load_system("test_tetrag_hfo2.xyz")

compute_local_descriptors(sys,ace)
compute_force_descriptors(sys,ace)
InteratomicPotentials.force(sys, lb)

f_check = Molly.forces(inter_ace, sys)

mp = molly_params(sys)
m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(typeof(1.0u"eV/Å"), 1),),
            force_units=u"eV/Å",
            energy_units=u"eV",
            #loggers=(force=ForceLogger(Float32, 1),),
            #force_units=NoUnits,
            )

simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(m_sys,simulator,0)
f0 = m_sys.loggers.force.history[1] # These forces match

simulate!(m_sys,simulator,100)

fcheck_log = [ustrip.(f) for f in m_sys.loggers.force.history[end]]
fcheck_ip = InteratomicPotentials.force(m_sys,lb)

atoms_arr = [AtomsBase.Atom(sp,pos) for (sp,pos) in zip(atomic_symbol(m_sys),position(m_sys))]
fs_check = FlexibleSystem(atoms_arr,bounding_box(m_sys),boundary_conditions(m_sys))

fcheck_fs = InteratomicPotentials.force(m_sys,lb)

function fcomp_errors(f_ref,f_test)
    comp_errs = []
    for i in 1:length(f_ref)
        for j in 1:3
            err = f_test[i][j]-f_ref[i][j]
            push!(comp_errs,err)
        end
    end
    comp_errs
end

@show maximum(fcomp_errors(fcheck_log,fcheck_ip)) # 0.0
@show maximum(fcomp_errors(fcheck_log,fcheck_fs)) # 0.0