#import Molly: forces, potential_energy
using Molly
using InteratomicPotentials
using Unitful, UnitfulAtomic
#import AtomsBase: AbstractSystem
using AtomsBase


e_dim = dimension(u"eV")
l_dim = dimension(u"Å")

def_eunit = u"eV"
def_lunit = u"Å"

nounit_t = typeof(NoUnits)

struct InteratomicPotentialInter{P<:AbstractPotential}
    potential::P
    energy_units::Unitful.Unitlike 
    length_units::Unitful.Unitlike

    # internal constructor, ensuring energy units and length units have correct dimensions
    ( InteratomicPotentialInter(pot::AbstractPotential, 
                                eu::Union{nounit_t, Unitful.Units{UE,e_dim,nothing}} = def_eunit, 
                                lu::Union{nounit_t, Unitful.Units{UL,l_dim,nothing}} = def_lunit) 
                                where {UE,UL} = new{typeof(pot)}(pot,eu,lu) )
end


function Molly.forces(inter::InteratomicPotentialInter, 
                      sys::AbstractSystem,
                      neighbors = nothing;
                      n_threads = Threads.nthreads())

    forces = InteratomicPotentials.force(sys,inter.potential)

    # initial profiling didn't show huge performance hit from unit conversion 
    # but! that may be because other parts of IP.jl are very slow, e.g. neighbor list construction
    if eltype(forces[1]) <: Unitful.Quantity 
        if inter.energy_units != NoUnits
            forces = [uconvert.(inter.energy_units/inter.length_units, fi)
                     for fi in forces]
        else
            forces = [ustrip.(fi)
                     for fi in forces]
        end
    elseif eltype(forces[1]) <: Real && inter.energy_units != NoUnits
        forces = forces * inter.energy_units/inter.length_units
    end

    forces
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

