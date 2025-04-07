#=
Expect this to move to IP.jl eventually
TODO:
-[ ] once everything satisfies updated AtomsCalcualtors interface, should be able to design constructor so that it autodetects member units and propagates them to CommitteePotential
=#
#export CommitteePotential

using InteratomicPotentials: AbstractPotential

const E_DIM = dimension(u"eV")
const EN_DIM = dimension(u"kJ * mol^-1") # need to account for annoying per atom case
const L_DIM = dimension(u"Å")
const NOUNIT_T = typeof(NoUnits)

#=
 - Assumes member potentials are same type
 - Because units are set indepently from members, risk of unexpected behavior where member potential units do not match committee potential units
=#
struct CommitteePotential{P}
    members::Vector{P}
    leader::Integer
    energy_units::Unitful.Unitlike
    length_units::Unitful.Unitlike
    force_units::Unitful.Unitlike

    #inner constructor to enforce that leader can appropriately access members
    function CommitteePotential(
                members::Vector{P},
                leader::Integer,
                eunits::Union{NOUNIT_T, Unitful.Units{UE,E_DIM,nothing},Unitful.Units{UE,EN_DIM,nothing}},
                lunits::Union{NOUNIT_T, Unitful.Units{UL,L_DIM,nothing}}) where {P, UE, UL}

        if !checkbounds(Bool, members, leader)
            error("leader index is out of bounds with respect to members array")
        end
        funits = eunits/lunits
        return new{P}(members,leader,eunits,lunits,funits)
    end
end

function CommitteePotential(members::Vector{P},
                            leader::Integer=1;
                            energy_units=u"eV",
                            length_units=u"Å") where {P}
    cmte_pot = CommitteePotential(members,leader,energy_units,length_units)
    cmte_pot
end

function AtomsCalculators.potential_energy(sys::AbstractSystem,
                                           cmte_pot::CommitteePotential;
                                           neighbors=nothing,
                                           n_threads::Integer=Threads.nthreads())
  # This manual type stuff is only needed so long as IP.jl doesn't satistfy the AtomsCalculators interface
  leader_pot_type = eltype(cmte_pot.members)
  if leader_pot_type <: AbstractPotential
    poteng = InteratomicPotentials.potential_energy(sys,cmte_pot.members[cmte_pot.leader]) * cmte_pot.energy_units
  else
    poteng = AtomsCalculators.potential_energy(sys,cmte_pot.members[cmte_pot.leader])
  end

  poteng
end

function AtomsCalculators.forces(sys::AbstractSystem,
                                 cmte_pot::CommitteePotential;
                                 neighbors=nothing,
                                 n_threads::Integer=Threads.nthreads())
  # This manual type stuff is only needed so long as IP.jl doesn't satistfy the AtomsCalculators interface
  leader_pot_type = eltype(cmte_pot.members)
  if leader_pot_type <: AbstractPotential
    force_units = cmte_pot.energy_units/cmte_pot.length_units
    forces = InteratomicPotentials.force(sys,cmte_pot.members[cmte_pot.leader]) * force_units
  else
    forces = AtomsCalculators.forces(sys,cmte_pot.members[cmte_pot.leader])
  end

  forces
end

function compute_all_energies(sys::AbstractSystem, cmte_pot::CommitteePotential)
  # This manual type stuff is only needed so long as IP.jl doesn't satistfy the AtomsCalculators interface

  leader_pot_type = eltype(cmte_pot.members)
  if leader_pot_type <: AbstractPotential
    all_potengs = [ InteratomicPotentials.potential_energy(sys,pot) * cmte_pot.energy_units
                    for pot in cmte_pot.members]
  else
    all_potengs = [AtomsCalculators.potential_energy(sys,pot)
                    for pot in cmte_pot.members]
  end

  all_potengs
end

function compute_all_forces(sys::AbstractSystem, cmte_pot::CommitteePotential)
  # This manual type stuff is only needed so long as IP.jl doesn't satistfy the AtomsCalculators interface

  leader_pot_type = eltype(cmte_pot.members)
  if leader_pot_type <: AbstractPotential
    force_units = cmte_pot.energy_units/cmte_pot.length_units
    all_forces = [ InteratomicPotentials.force(sys,pot) * force_units
                  for pot in cmte_pot.members]
  else
    all_forces = [ AtomsCalculators.forces(sys,pot)
                  for pot in cmte_pot.members]
  end

  all_forces
end


function get_params(cmte_pot::CommitteePotential)
  all_params = Float64[]
  for member in cmte_pot.members
    params = member.params
    push!(all_params,params)
  end
end
