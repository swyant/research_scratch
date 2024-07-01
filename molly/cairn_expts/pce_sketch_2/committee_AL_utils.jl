import InteratomicPotentials: potential_energy, force
using AtomsBase
#=
TODO Will need to make this AtomsCalculators compatible, but I need to do that in IP.jl first
=#

struct CommitteePotential 
  members::Vector{<:AbstractPotential} #Should be AbstractPotential
  leader::Integer
end 

function CommitteePotential(members::Vector{<:AbstractPotential})
    cmte_pot = CommitteePotential(members,1)
    cmte_pot
end

function potential_energy(sys::AbstractSystem, cmte_pot::CommitteePotential)
  poteng = potential_energy(sys,cmte_pot.members[cmte_pot.leader])
  poteng
end

function force(sys::AbstractSystem, cmte_pot::CommitteePotential)
  forces = force(sys,cmte_pot.members[cmte_pot.leader])
  forces
end

function compute_all_energies(sys::AbstractSystem, cmte_pot::CommitteePotential)
  all_potengs = [ potential_energy(sys,pot)
                  for pot in cmte_pot.members]
  all_potengs
end

function compute_all_forces(sys::AbstractSystem, cmte_pot::CommitteePotential)
  all_forces = [ force(sys,pot)
                for pot in cmte_pot.members]
  all_forces
end


##-----##

abstract type AbstractCommitteeQoI end

struct CommitteeEnergy <: AbstractCommitteeQoI
    cmte_reduce::Union{Nothing,Function}
end

function CommitteeEnergy()
  cmte_e_qoi = CommitteeEnergy(nothing)
  cmte_e_qoi
end

function compute(qoi::CommitteeEnergy,sys::AbstractSystem,cmte_pot::CommitteePotential)
    energies = compute_all_energies(sys,cmte_pot)
    if !isnothing(qoi.cmte_reduce)
      reduced_energies = qoi.cmte_reduce(energies)
      return reduced_energies
    else
      return energies
    end
end

#TODO need to be able to return site energies from arbitrary AbstractPotential
#TODO inner constructor to enforce no (nothing,something)
#struct CommitteeAtomicEnergy <: AbstractCommitteeQoI
#    cmte_reduce::Union{Nothing,Function}    
#    atom_reduce::Union{Nothing,Function}
#end

#TODO subselect number of atoms
struct CommitteeFlatForces <: AbstractCommitteeQoI
    cmte_reduce::Union{Nothing,Function}
    coord_and_atom_reduce::Union{Nothing,Function}
    reduce_order::Vector{Int64}
end

# only run @ initialization of committee QoI
function _check_reduction_fn(fn::Function)
  test_arr1 = [1.0,2.0,3.0,4.0,5.0]

  local res = nothing # to deal with the scoping in the try statement, can now use res in else
  try
    res = fn(test_arr1)
  catch 
  else
    check_res = typeof(res) <: Union{<:Real, <:Integer, Bool} 
    return check_res
  end

  test_arr2 = [1,2,3,4,5]
  try
    res = fn(test_arr2)
  catch 
  else
    check_res = typeof(res) <: Union{<:Real, <:Integer, Bool} 
    return check_res
  end

  test_arr3 = [true,true,false,true,false]
  try
    res = fn(test_arr3)
  catch 
  else
    check_res = typeof(res) <: Union{<:Real, <:Integer, Bool} 
    return check_res
  end

  false
end

function CommitteeFlatForces(nt::NamedTuple{<:Any, <:Tuple{Vararg{Function}}})
    fn_dict = Dict(:cmte => 
                    Dict{String,Union{Nothing,Function,Int64}}(
                        "fn"  => nothing,
                        "idx" => 1),
                    :coord_and_atom =>
                    Dict{String,Union{Nothing,Function,Int64}}(
                        "fn"  => nothing,
                        "idx" => 2)
                    )

    if !all(in.(keys(nt),[keys(fn_dict)]))
      error("""Only allowed keys are "cmte", "coord_and_atom" """)
    elseif length(nt) > 2
      error("There can be a maximum of 2 elements in the passed NamedTuple")
    end
    
    #TODO this should really be in the inner constructor
    for fn in nt
      if !_check_reduction_fn(fn)
        error("reduction function must reduce to single Float, Integer, or Boolean")
      end
    end 

    reduce_order = Int64[]
    for (k,fn) in pairs(nt)
      push!(reduce_order,fn_dict[k]["idx"])
      fn_dict[k]["fn"] = fn
    end

    cmte_flat_force_qoi = CommitteeFlatForces(fn_dict[:cmte]["fn"],
                                         fn_dict[:coord_and_atom]["fn"],
                                         reduce_order)
    cmte_flat_force_qoi
end

function CommitteeFlatForces()
    flat_force_qoi = CommitteeFlatForces(nothing,nothing,Int64[])
    flat_force_qoi
end

function compute(qoi::CommitteeFlatForces,sys::AbstractSystem,cmte_pot::CommitteePotential)
    reduce_dict = Dict{Int64, Union{Nothing,Function}}(
                  1 => qoi.cmte_reduce,
                  2 => qoi.coord_and_atom_reduce)


    raw_all_forces = compute_all_forces(sys,cmte_pot)
    all_forces = stack(map(elem->stack(elem,dims=1),raw_all_forces), dims=1)  # num_cmte x num_atoms x 3
    all_flat_forces = reshape(permutedims(all_forces,(1,3,2)), size(all_forces,1), :) #num_cmte x num_atoms*3 1x,1y,1z,2x,2y,2z,etc.

    if isnothing(qoi.cmte_reduce) && isnothing(qoi.coord_and_atom_reduce)
      return all_flat_forces
    else 
      inter = all_flat_forces 
      for d in qoi.reduce_order
        inter = mapslices(reduce_dict[d],inter,dims=d)
      end

      if length(qoi.reduce_order) == 2 
          @assert size(inter) == (1,1) && typeof(inter[1]) <: Union{<:Real, <:Integer, Bool}
          final_qoi = inter[1]
      else 
          for d in qoi.reduce_order
            @assert size(inter,d) == 1 # ensure singleton dimension
          end 
          
          #arguably should check if Int,bool, float, but once I put that in the inner constructor it's fine
          final_qoi = dropdims(inter,dims=Tuple(qoi.reduce_order))
      end 
      return final_qoi
    end
end


struct CommitteeForces <: AbstractCommitteeQoI
    cmte_reduce::Union{Nothing,Function}
    atom_reduce::Union{Nothing,Function}
    coord_reduce::Union{Nothing,Function}
    reduce_order::Vector{Int64}
end

function CommitteeForces()
  cmte_force_qoi = CommitteeForces(nothing,nothing,nothing,Int64[])
  cmte_force_qoi
end

function CommitteeForces(nt::NamedTuple{<:Any, <:Tuple{Vararg{Function}}})
    fn_dict = Dict(:cmte => 
                    Dict{String,Union{Nothing,Function,Int64}}(
                        "fn"  => nothing,
                        "idx" => 1),
                    :atom =>
                    Dict{String,Union{Nothing,Function,Int64}}(
                        "fn"  => nothing,
                        "idx" => 2),
                    :coord =>
                    Dict{String,Union{Nothing,Function,Int64}}(
                        "fn"  => nothing,
                        "idx" => 3),

                    )

    if !all(in.(keys(nt),[keys(fn_dict)]))
      error("""Only allowed keys are "cmte", "atom", "coord" """)
    elseif length(nt) > 3
      error("There can be a maximum of 3 elements in the passed NamedTuple")
    end
    
    #TODO this should really be in the inner constructor
    for fn in nt
      if !_check_reduction_fn(fn)
        error("reduction function must reduce to single Float, Integer, or Boolean")
      end
    end 

    reduce_order = Int64[]
    for (k,fn) in pairs(nt)
      push!(reduce_order,fn_dict[k]["idx"])
      fn_dict[k]["fn"] = fn
    end

    cmte_force_qoi = CommitteeForces( fn_dict[:cmte]["fn"],
                                      fn_dict[:atom]["fn"],
                                      fn_dict[:coord]["fn"],
                                      reduce_order)
    cmte_force_qoi
end

function compute(qoi::CommitteeForces,sys::AbstractSystem,cmte_pot::CommitteePotential)
  reduce_dict = Dict{Int64, Union{Nothing,Function}}(
                1 => qoi.cmte_reduce,
                2 => qoi.atom_reduce,
                3 => qoi.coord_reduce)

  raw_all_forces = compute_all_forces(sys,cmte_pot)

  if isnothing(qoi.cmte_reduce) && isnothing(qoi.atom_reduce) && isnothing(qoi.coord_reduce)
    return all_flat_forces
  else 
    inter = stack(map(elem->stack(elem,dims=1),raw_all_forces), dims=1)  # num_cmte x num_atoms x 3
    for d in qoi.reduce_order
      inter = mapslices(reduce_dict[d],inter,dims=d)
    end

    if length(qoi.reduce_order) == 3 
        @assert size(inter) == (1,1,1) && typeof(inter[1]) <: Union{<:Real, <:Integer, Bool}
        final_qoi = inter[1]
    else 
        for d in qoi.reduce_order
          @assert size(inter,d) == 1 # ensure singleton dimension
        end 
        
        #arguably should check if Int,bool, float, but once I put that in the inner constructor it's fine
        final_qoi = dropdims(inter,dims=Tuple(qoi.reduce_order))
    end 
    return final_qoi
  end
end