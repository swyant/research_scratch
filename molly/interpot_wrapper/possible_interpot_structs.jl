using InteratomicPotentials
using Unitful, UnitfulAtomic
#=
Goals: 
1. Want energy units and force units (or length units) to be appropriately type constrained 
2. Ideally, InteratomicPotentialInter is only parameterized on the abstractpotential 
=#

e_dim = dimension(u"eV")
f_dim = dimension(u"eV/Å")


# This is OK, but I'd prefer not to have to parameterize the type off of UE, UF as well. 
#struct InteratomicPotentialInter{P<:AbstractPotential, UE, UF}
#    potential::P
#    energy_units::Unitful.Units{UE,e_dim,nothing}
#    force_units::Unitful.Units{UF,f_dim,nothing}
#end

# similar to https://github.com/ACEsuit/ACEmd.jl/blob/main/src/structs.jl 
#struct InteratomicPotentialInter{P<:AbstractPotential}
#    potential::P
#    energy_units::Unitful.Unitlike
#    force_units::Unitful.Unitlike # don't think it's actually necessarry to type these fields if relying on inner constructor?
#
#    function InteratomicPotentialInter(pot::AbstractPotential, 
#                                       eu::Unitful.Unitlike, 
#                                       fu::Unitful.Unitlike)
#        @assert dimension(eu) == e_dim 
#        @assert dimension(fu) == f_dim
#        new{typeof(pot)}(pot,eu,fu)
#    end
#end



# my preferred solution
#struct InteratomicPotentialInter{P<:AbstractPotential}
#    potential::P
#    energy_units::Unitful.Unitlike
#    force_units::Unitful.Unitlike
#
#    InteratomicPotentialInter(pot::AbstractPotential, 
#                              eu::Unitful.Units{UE,e_dim,nothing}, 
#                              fu::Unitful.Units{UF,f_dim,nothing}) where {UE,UF} = new{typeof(pot)}(pot,eu,fu)
#end


## Should I also enforce defaults in the inner constructor like they do in AceMD?
def_eunit = u"eV"
def_funit = u"eV/Å"

struct InteratomicPotentialInter{P<:AbstractPotential}
    potential::P
    energy_units::Unitful.Unitlike
    force_units::Unitful.Unitlike

    InteratomicPotentialInter(pot::AbstractPotential, 
                              eu::Unitful.Units{UE,e_dim,nothing} = def_eunit, 
                              fu::Unitful.Units{UF,f_dim,nothing} = def_funit) where {UE,UF} = new{typeof(pot)}(pot,eu,fu)
end
