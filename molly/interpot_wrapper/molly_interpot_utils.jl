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

