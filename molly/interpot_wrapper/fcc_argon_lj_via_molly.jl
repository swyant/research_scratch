using Molly
using InteratomicPotentials
using AtomsIO 
using Unitful
using AtomsBase
 
# Wrapper struct enabling potentials implemented in InteratomicPotentials.jl to be used in Molly
struct InteratomicPotentialInter{P<:AbstractPotential}
	potential::P
    eaf::Function
    # TODO: should also use the energy cache Jeremiah had, but ignoring that for now
    # TODO: Default constructur with energy_and_force probably
end

# Extending Molly.forces with above struct
function Molly.forces(inter::InteratomicPotentialInter, sys::AbstractSystem, neighbors=nothing; n_threads=Threads.nthreads())
    eandf = inter.eaf(sys, inter.potential)
    eandf.f
end

function molly_params(sys::AtomsBase.AbstractSystem)
    coords = [SVector{3}(pos) for pos in position(sys)] # need to be SVector for zero() func to work
    atoms = [Molly.Atom(mass=atm_mass) for atm_mass in atomic_mass(sys)]
    atoms_data = [AtomData(element=string(atm_symbol)) for atm_symbol in atomic_symbol(sys)]
    velocities = [SVector{3}(vel) for vel in velocity(sys)] 

    # TODO: Generalize this. Currently very fragile, assumes you read in a cubic box. Otherwise need to use TriclinicBoundary, but some limitations...
    # Also by default, assumes directions are periodic..., ignoring boundary_conditions
    bbox = bounding_box(sys)
    boundary=CubicBoundary(bbox[1][1],bbox[2][2],bbox[3][3]) # FRAGILE! 

    return NamedTuple((atoms=atoms, atoms_data=atoms_data, coords=coords, velocities = velocities, boundary=boundary))
end 

############################################################

# load system with AtomsIO
sys = load_system("../ref_config/dump_final.extxyz")


# set up InteratomicPotentials LennardJones, based off of 10.1103/PhysRevB.54.340 
ϵ = 0.01032u"eV"
σ = 3.405u"Å"
rcut = 8.51u"Å"
species = [:Ar]
lj_p = InteratomicPotentials.LennardJones(ϵ,σ,rcut,species)

inter_lj = InteratomicPotentialInter(lj_p,InteratomicPotentials.energy_and_force)
general_inters = (inter_lj,)

# Check force with regular sys obtained with AtomsIO
f_p = Molly.forces(inter_lj, sys)
#f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]

mp = molly_params(sys)

m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(typeof(1.0u"eV/Å"), 1),),
            force_units=u"eV/Å",
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
lines = readlines("../ref_config/raw_forces")
f_ref = [parse.(Float64,split(li)) for li in lines]u"eV/Å"

fcomp_errs = []
for i in 1:length(f_ref)
    for j in 1:3
        err = f_check[i][j]-f_ref[i][j]
        push!(fcomp_errs,err)
    end
end

@show maximum(fcomp_errs)