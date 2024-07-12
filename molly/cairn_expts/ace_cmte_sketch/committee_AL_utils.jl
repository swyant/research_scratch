import InteratomicPotentials: InteratomicPotentials, AbstractPotential
import AtomsCalculators
using AtomsBase
import Molly: Molly, System, log_property!  # do I need the Molly?
using Unitful
using Cairn

#=
TODO Will need to make this (properly!)AtomsCalculators compatible, but I need to do that in IP.jl first
Right now it's just hacky
=#

struct CommitteePotential{F,E}
  members::Union{Vector{<:AbstractPotential}, Vector{<:PolynomialChaos}}#Should be AbstractPotential
  leader::Int64
  force_units::F 
  energy_units::E
end 

function CommitteePotential(members::Union{Vector{<:AbstractPotential}, Vector{<:PolynomialChaos}},
                            leader=1;
                            force_units=u"eV/Å",
                            energy_units=u"eV")
    cmte_pot = CommitteePotential(members,leader,force_units,energy_units)
    cmte_pot
end

function AtomsCalculators.potential_energy(sys::AbstractSystem, 
                                           cmte_pot::CommitteePotential;
                                           neighbors=nothing,
                                           n_threads::Integer=Threads.nthreads())
  # This manual type stuff is only needed so long as IP.jl doesn't satistfy the AtomsCalculators interface
  leader_pot_type = typeof(cmte_pot.members[cmte_pot.leader])
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
  leader_pot_type = typeof(cmte_pot.members[cmte_pot.leader])
  if leader_pot_type <: AbstractPotential
    forces = InteratomicPotentials.force(sys,cmte_pot.members[cmte_pot.leader]) * cmte_pot.force_units
  else 
    forces = AtomsCalculators.forces(sys,cmte_pot.members[cmte_pot.leader])
  end

  forces
end

function compute_all_energies(sys::AbstractSystem, cmte_pot::CommitteePotential)
  # This manual type stuff is only needed so long as IP.jl doesn't satistfy the AtomsCalculators interface
  leader_pot_type = typeof(cmte_pot.members[cmte_pot.leader])
  if leader_pot_type <: AbstractPotential
    all_potengs = [ InteratomicPotentials.potential_energy(sys,pot)
                    for pot in cmte_pot.members]
  else
    all_potengs = [AtomsCalculators.potential_energy(sys,pot)
                    for pot in cmte_pot.members]
  end
  all_potengs
end

function compute_all_forces(sys::AbstractSystem, cmte_pot::CommitteePotential)
  # This manual type stuff is only needed so long as IP.jl doesn't satistfy the AtomsCalculators interface
  leader_pot_type = typeof(cmte_pot.members[cmte_pot.leader])
  if leader_pot_type <: AbstractPotential
    all_forces = [ InteratomicPotentials.force(sys,pot)
                  for pot in cmte_pot.members]
  else 
    all_forces = [ AtomsCalculators.forces(sys,pot)
                  for pot in cmte_pot.members]
  end                
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

function compute(qoi::CommitteeEnergy,sys::AbstractSystem,cmte_pot::CommitteePotential; cache_field=nothing)
    
    if typeof(sys) <: Molly.System && !isnothing(cache_field)
      if isnothing(sys.data[cache_field])
        all_energies = compute_all_energies(sys,cmte_pot)
        sys.data[cache_field] = all_energies 
      else
        all_energies = sys.data[cache_field]
      end
    else
        all_energies = compute_all_energies(sys,cmte_pot)
    end

    if !isnothing(qoi.cmte_reduce)
      reduced_energies = qoi.cmte_reduce(all_energies)
      return reduced_energies
    else
      return all_energies
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

function compute(qoi::CommitteeFlatForces,sys::AbstractSystem,cmte_pot::CommitteePotential; cache_field=nothing)
    reduce_dict = Dict{Int64, Union{Nothing,Function}}(
                  1 => qoi.cmte_reduce,
                  2 => qoi.coord_and_atom_reduce)

    if typeof(sys) <: Molly.System && !isnothing(cache_field)
      if isnothing(sys.data[cache_field])
        raw_all_forces = compute_all_forces(sys,cmte_pot)
        sys.data[cache_field] = raw_all_forces
      else
        raw_all_forces = sys.data[cache_field]
      end
    else
        raw_all_forces = compute_all_forces(sys,cmte_pot)
    end

    #TODO Is it worth it to cache flat forces directly rather than raw_all_forces
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
        # assert statement fails with units 
        #@assert size(inter) == (1,1) && typeof(inter[1]) <: Union{<:Real, <:Integer, Bool}
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

function compute(qoi::CommitteeForces,sys::AbstractSystem,cmte_pot::CommitteePotential; cache_field=nothing)
  reduce_dict = Dict{Int64, Union{Nothing,Function}}(
                1 => qoi.cmte_reduce,
                2 => qoi.atom_reduce,
                3 => qoi.coord_reduce)

  if typeof(sys) <: Molly.System && !isnothing(cache_field)
    if isnothing(sys.data[cache_field])
      raw_all_forces = compute_all_forces(sys,cmte_pot)
      sys.data[cache_field] = raw_all_forces
    else
      raw_all_forces = sys.data[cache_field]
    end
  else
      raw_all_forces = compute_all_forces(sys,cmte_pot)
  end

  if isnothing(qoi.cmte_reduce) && isnothing(qoi.atom_reduce) && isnothing(qoi.coord_reduce)
    return raw_all_forces
  else 
    inter = stack(map(elem->stack(elem,dims=1),raw_all_forces), dims=1)  # num_cmte x num_atoms x 3
    for d in qoi.reduce_order
      inter = mapslices(reduce_dict[d],inter,dims=d)
    end

    if length(qoi.reduce_order) == 3 
        # assert statement fails with units 
        #@assert size(inter) == (1,1,1) && typeof(inter[1]) <: Union{<:Real, <:Integer, Bool}
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


##========###
mutable struct SimpleTriggerLogger{T}
  observable::Union{Nothing,T}
  n_steps::Int
  history::Vector{T}
end

# if n_steps=0, only observable is recorded for the length of that timestep
function SimpleTriggerLogger(n_steps::Integer=0, T::DataType=Float64) 
  return SimpleTriggerLogger{T}(nothing, n_steps, T[])
end

Base.values(logger::SimpleTriggerLogger) = logger.history

# Basically skip the standard log_property!() call
function log_property!(logger::SimpleTriggerLogger, s::System, neighbors=nothing,
  step_n::Integer=0; n_threads::Integer=Threads.nthreads(), kwargs...)
end

function log_property!(logger::SimpleTriggerLogger{T}, obs::T, step_n::Integer=0) where {T}

  logger.observable = obs

  if logger.n_steps != 0 && (step_n % logger.n_steps) == 0
        push!(logger.history, obs)
  end
end

function reset_observable!(logger::SimpleTriggerLogger)
  logger.observable = nothing
end

#function Base.show(io::IO, fl::SimpleTriggerLogger)
#  print(io, "SimpleTriggerLogger{}", eltype(fl.trigger), ", ", eltype(eltype(values(fl))), "} with n_steps ",
#          fl.n_steps, ", ", length(values(fl)), " frames recorded for ",
#          length(values(fl)) > 0 ? length(first(values(fl))) : "?", " atoms")
#end


##========###

function initialize_triggers(triggers::Tuple, sys::Molly.System) 
  if typeof(sys.loggers) <: Tuple && length(sys.loggers) == 0  # if loggers empty Tuple, convert to empty NamedTuple
    loggers = NamedTuple()
  elseif typeof(sys.loggers) <: NamedTuple
    loggers = sys.loggers
  else
    error("oops, I can't handle a case where sys.loggers is Tuple with finite size, or anything other than NamedTuple")
  end

  if isnothing(sys.data)
    ddict = Dict{Any,Any}() #have to be flexible with types, user can do anything
  elseif sys.data <: Dict
    ddict = Dict{Any,Any}(sys.data) # existing data dict may be too strictly typed
  else 
    error("System.data needs to be either nothing or a dictionary")
  end

  for trigger in triggers
    loggers = append_loggers(trigger,loggers)
    ddict   = initialize_data(trigger,ddict)
  end

  return_sys = Molly.System(sys; loggers=loggers, data=ddict)
  #return_sys = Molly.System(sys; loggers=loggers)

  return_sys
end

struct CmteTrigger{T <:AbstractCommitteeQoI, F<:Function} <: ActiveLearningTrigger
  cmte_qoi::T
  compare::F
  thresh::Union{Float64,Int64,Bool} #worth it to generalize?
  cmte_pot::Union{Nothing,CommitteePotential}
  logger_spec::Union{Nothing,Tuple{Symbol,Int64}}
end 

function CmteTrigger(cmte_qoi::AbstractCommitteeQoI,
                     compare::Function, 
                     thresh::Union{Float64,Int64,Bool},
                     cmte_pot::Union{Nothing,CommitteePotential}=nothing;
                     logger_spec::Union{Nothing,Tuple{Symbol,Int64}}=nothing)
  cmte_trigger = CmteTrigger(cmte_qoi,compare,thresh,cmte_pot,logger_spec)
  cmte_trigger
end

# Should this modify the system directly, return a modified loggers tuple, or return something like (new_keys), (new_loggers)
function append_loggers(trigger::CmteTrigger, loggers::NamedTuple)

  if !isnothing(trigger.logger_spec)
    prior_keys = keys(loggers)
    prior_vals = values(loggers) 

    log_symb = trigger.logger_spec[1]
    @assert !(log_symb in prior_keys)
    log_freq = trigger.logger_spec[2]
    log_type = typeof(trigger.thresh) #TODO: need a more robust way of typing the CommitteeQoIs

    new_logger = SimpleTriggerLogger(log_freq,log_type)

    updated_keys = (prior_keys...,log_symb)
    updated_vals = (prior_vals...,new_logger)

    return NamedTuple{updated_keys}(updated_vals)
  else 
    return loggers
  end
end

function initialize_data(trigger::CmteTrigger, ddict::Dict)
  return ddict
end 

function trigger_activated!(trigger::CmteTrigger, 
                            sys::Molly.System, 
                            step_n::Integer=0; 
                            shared_cmte_pot::Union{Nothing,CommitteePotential}=nothing,
                            cache_field=nothing)

    #Select commitee potential : shared_cmte_pot > trigger.cmte_pot > sys.general_inters[1]
    if !isnothing(shared_cmte_pot)
      cmte_pot = shared_cmte_pot # prefer over CmteTrigger.cmte_pot
    elseif !isnothing(trigger.cmte_pot)
      cmte_pot = trigger.cmte_pot
    elseif typeof(sys.general_inters[1]) <: CommitteePotential
      cmte_pot = sys.general_inters[1]
    else 
      error("No committee potential available for trigger activation")
    end

    # Assuming it has to be unitless right now
    cmte_qoi = ustrip(compute(trigger.cmte_qoi,sys,cmte_pot;cache_field=cache_field))
    @show cmte_qoi
    @assert typeof(ustrip(cmte_qoi)) <: typeof(trigger.thresh)
    
    if !isnothing(trigger.logger_spec)
      log_property!(sys.loggers[trigger.logger_spec[1]], cmte_qoi,step_n)
    end

    res = trigger.compare(cmte_qoi,trigger.thresh)
    res
end

struct SharedCmteTrigger <: ActiveLearningTrigger
  cmte_pot::CommitteePotential # technically this can be nothing if user wants to use sys.general_inters[1]
  subtriggers::Tuple{Vararg{<:CmteTrigger}} # Fundamentally a bit type unstable
  energy_cache_field::Union{Nothing,Symbol}
  force_cache_field::Union{Nothing,Symbol}
end

#TODO this should just be inner constructor maybe?
function SharedCmteTrigger(cmte_pot::CommitteePotential, 
                           subtriggers::Tuple{Vararg{<:CmteTrigger}};
                           energy_cache_field = nothing,
                           force_cache_field = nothing)
  shared_cmte_trigger = SharedCmteTrigger(cmte_pot,subtriggers,energy_cache_field, force_cache_field)
  shared_cmte_trigger
end

function append_loggers(shared_trigger::SharedCmteTrigger, loggers::NamedTuple)
  for subtrigger in shared_trigger.subtriggers
    loggers = append_loggers(subtrigger,loggers)
  end
  return loggers
end


function initialize_data(shared_trigger::SharedCmteTrigger, ddict::Dict)
  ecache_field = shared_trigger.energy_cache_field
  fcache_field = shared_trigger.force_cache_field

  if !isnothing(ecache_field) || !isnothing(fcache_field)
    if !(:_reset_every_step in keys(ddict))
        ddict[:_reset_every_step] = Symbol[]
    end

    if !isnothing(ecache_field)
      @assert !(ecache_field in keys(ddict))
      ddict[ecache_field] = nothing
      push!(ddict[:_reset_every_step],ecache_field)
    end

    if !isnothing(fcache_field)
      @assert !(fcache_field in keys(ddict))
      ddict[fcache_field] = nothing
      push!(ddict[:_reset_every_step],fcache_field)
    end

    return ddict
  else
    return ddict
  end
end 

# not handling the case where the shared cmte pot is 
function trigger_activated!(shared_trigger::SharedCmteTrigger,sys::Molly.System, step_n::Integer=0)
  all_res = Bool[]
  for subtrigger in shared_trigger.subtriggers
    #How to pass appropriate cache field for custom Committee QoIs
    if typeof(subtrigger.cmte_qoi) <: Union{CommitteeForces, CommitteeFlatForces}
      res = trigger_activated!(subtrigger,sys, step_n; shared_cmte_pot=shared_trigger.cmte_pot, cache_field=shared_trigger.force_cache_field) 
    elseif typeof(subtrigger.cmte_qoi) <: CommitteeEnergy
      res = trigger_activated!(subtrigger,sys, step_n; shared_cmte_pot=shared_trigger.cmte_pot, cache_field=shared_trigger.energy_cache_field) 
    else 
      res = trigger_activated!(subtrigger,sys, step_n; shared_cmte_pot=shared_trigger.cmte_pot) 
    end
    push!(all_res, res)
  end

  final_res = any(all_res)
  final_res
end

#=
This should be a reset_per_step!(sys, triggers) where you loop over triggers, and have a new reset_trigger_observables()
=#

function reset_logger!(trigger::CmteTrigger, sys::Molly.System)
  if !isnothing(trigger.logger_spec)
    reset_observable!(sys.loggers[trigger.logger_spec[1]])
  end
end

function reset_logger!(shared_trigger::SharedCmteTrigger, sys::Molly.System)
  for subtrigger in shared_trigger.subtriggers
    reset_logger!(subtrigger,sys)
  end
end

function perstep_reset!(triggers, sys::Molly.System)
  # need an API to get whether trigger has associated logger and what the logger name is
  for trigger in triggers 
    reset_logger!(trigger,sys)
  end

  if (typeof(sys.data) <: Dict && 
    :_reset_every_step in keys(sys.data) && 
    length(sys.data[:_reset_every_step]) > 0)

    for dict_symb in sys.data[:_reset_every_step] 
      sys.data[dict_symb] = nothing 
    end
  end
end
