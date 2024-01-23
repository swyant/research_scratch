using AtomsIO

include("./alt_interpot_utils.jl")

# load system with AtomsIO
sys = load_system(ExtxyzParser(), "dump_final.xyz")


# set up InteratomicPotentials LennardJones, based off of 10.1103/PhysRevB.54.340 
ϵ = 0.01032u"eV"
σ = 3.405u"Å"
rcut = 8.51u"Å"
species = [:Ar]
lj_p = InteratomicPotentials.LennardJones(ϵ,σ,rcut,species)

general_inters = (lj_p,)

# Check force with regular sys obtained with AtomsIO
f_p = Molly.forces(lj_p, sys)
f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]

mp = molly_params(sys)

m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(typeof(1.0u"eV/Å"), 1),),
            force_units=u"eV/Å",
            energy_units=u"eV",
            #loggers=(force=ForceLogger(Float32, 1),),
            #force_units=NoUnits,
            )

#### run zero simulation
simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(m_sys,simulator,0)

f_check = [f for f in m_sys.loggers.force.history[1]]

# Check the forces 
lines = readlines("./ref_forces")
f_ref = [parse.(Float64,split(li)) for li in lines]u"eV/Å"

fcomp_errs = []
for i in 1:length(f_ref)
    for j in 1:3
        err = f_check[i][j]-f_ref[i][j]
        push!(fcomp_errs,err)
    end
end

@show maximum(fcomp_errs)
# 7.563394355258879e-15 eV Å⁻¹
