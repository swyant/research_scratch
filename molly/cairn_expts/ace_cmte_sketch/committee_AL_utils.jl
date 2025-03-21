import InteratomicPotentials: InteratomicPotentials, AbstractPotential
import AtomsCalculators
using AtomsBase
import Molly: Molly, System, log_property!  # do I need the Molly?
import PotentialLearning: SubsetSelector, AbstractLearningProblem, learn!
using Unitful
using Cairn

#=
TODO Will need to make this (properly!)AtomsCalculators compatible, but I need to do that in IP.jl first
Right now it's just hacky
=#

# maybe reorder energy and force units
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


function get_params(mlip)
  params = mlip.params
end

function get_params(cmte_pot::CommitteePotential)
  all_params = Float64[]
  for member in cmte_pot.members
    params = member.params
    push!(all_params,params)
  end
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

# need to use this trick https://stackoverflow.com/questions/40160120/generic-constructors-for-subtypes-of-an-abstract-type
function (::Type{CmteTrigger{T, F}})(trigger::CmteTrigger{T,F};
                     cmte_qoi::T=trigger.cmte_qoi,
                     compare::F=trigger.compare,
                     thresh::Union{Float64,Int64,Bool}=trigger.thresh,
                     cmte_pot::Union{Nothing,CommitteePotential}=trigger.cmte_pot,
                     logger_spec::Union{Nothing,Tuple{Symbol,Int64}}=trigger.logger_spec) where {T<:AbstractCommitteeQoI, F<:Function}
  cmte_trigger = CmteTrigger{T,F}(trigger.cmte_qoi,trigger.compare,trigger.thresh,cmte_pot,logger_spec)
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

# also added the al argument, not sure if it should be there
function trigger_activated!(trigger::CmteTrigger,
                            sys::Molly.System,
                            al,
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


function SharedCmteTrigger(trigger::SharedCmteTrigger;
                           cmte_pot::CommitteePotential=trigger.cmte_pot,
                           subtriggers::Tuple{Vararg{<:CmteTrigger}}=trigger.subtriggers,
                           energy_cache_field = trigger.energy_cache_field,
                           force_cache_field = trigger.force_cache_field)
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

# not handling the case where the shared cmte pot is (what? spencer what? lol)
# also added the al argument, not sure if it should be there
function trigger_activated!(shared_trigger::SharedCmteTrigger,sys::Molly.System, al, step_n::Integer=0)
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

function get_logger_ids(shared_trigger::SharedCmteTrigger)
  trigger_ids = [get_logger_ids(subtrigger;from_shared=true)
                for subtrigger in shared_trigger.subtriggers]
  Tuple(trigger_ids)
end

function get_logger_ids(trigger::CmteTrigger; from_shared=false)
  if !isnothing(trigger.logger_spec)
    if from_shared
      return trigger.logger_spec[1]
    else
      return (trigger.logger_spec[1],)
    end
  else
    return nothing
  end
end

# Joanna's current way, immediate return true once any is satisfied
function trigger_activated(triggers::Tuple{Vararg{<:ActiveLearningTrigger}}, sys::Molly.System, al, step_n::Integer=1)
  for trigger in triggers
    if trigger_activated!(trigger, sys, al, step_n)
        return true
    end
end
end


####==========#####

mutable struct ALRoutine
  ref
  mlip
  trainset::Vector{<:AbstractSystem}
  triggers::Tuple{Vararg{<:ActiveLearningTrigger}}
  ss::SubsetSelector
  lp::AbstractLearningProblem #need to enforce concrete
  #trigger_updates::Union{Nothing,<:ALTriggerUpdate}
  trigger_updates
  aldata_spec
  cache::Dict
end

function initialize_al_cache!(al::ALRoutine)
  al.cache[:step_n] = nothing
  al.cache[:trainset_changes] = nothing
end

function strip_system(sys::Molly.System)
  stripped_sys = System(
                  atoms    = sys.atoms,
                  coords   = sys.coords,
                  boundary = sys.boundary)
  stripped_sys
end

#=
Should have an abstraction where for any kind of selector where you just append a new (labeled) configurated
=#

struct GreedySelector <: SubsetSelector
end

function update_trainset!(ss::GreedySelector,
                          sys::Molly.System,
                          al::ALRoutine)
  new_sys = strip_system(sys)
  new_trainset = reduce(vcat, [al.trainset, new_sys])

  new_trainset, [new_sys,]
end

# Recomputing all descriptors, energies, forces for entire trainset
struct InefficientLearningProblem <: AbstractLearningProblem
    weights::Vector{Float64}
    intcpt::Bool
    ref
end

function InefficientLearningProblem(weights=[1000.0,1.0],intcpt=false; ref=nothing)
    return InefficientLearningProblem(weights,intcpt,ref)
end

function learn(ilp::InefficientLearningProblem, mlip, trainset; ref=ilp.ref)
  lp = learn!(trainset, ref, mlip, ilp.weights, ilp.intcpt; e_flag=true, f_flag=true)
  new_mlip = deepcopy(mlip) #How to generalize? Should mlip be modified in place
  new_mlip.params = lp.β

  new_mlip
end

function retrain!(ilp::InefficientLearningProblem, sys::Molly.System, al::ALRoutine)
  # if this is a common pattern, then could just change the method handle in the AL loop
  # or have a generalized retrain!() that calls a more detailed retrain!() function
  trainset = al.trainset
  ref_pot = al.ref
  mlip = al.mlip

  new_mlip = learn(ilp,mlip,trainset;ref=ref_pot)

  new_mlip
end

####=====####

abstract type AbstractCmteLearningProblem <: AbstractLearningProblem end

# assumming only a simple append mode here. Original name: LinearCmteTriggerUpdate_Append
# Honestly there should be a cmtelearningproblem
mutable struct SubsampleAppendCmteRetrain  <: AbstractCmteLearningProblem
  lp::AbstractLearningProblem #need to enforce concrete
  cmte_indices::Vector{Vector{Integer}}
  #trainset (this could be a field, somewhat justified given tight coupling between indices and corresponding dataset)
end

# nothing is modified in place here
function learn(clp::SubsampleAppendCmteRetrain, old_cmte_pot, new_trainset)
  member_type = eltype(old_cmte_pot.members)
  new_members = Vector{member_type}()
  for (old_mlip, train_indices) in zip(old_cmte_pot.members, clp.cmte_indices)
    new_mlip = learn(clp.lp, old_mlip, new_trainset[train_indices])
    push!(new_members,new_mlip)
  end

  new_cmte_pot = CommitteePotential(new_members,
                                    old_cmte_pot.leader,
                                    old_cmte_pot.force_units,
                                    old_cmte_pot.energy_units)
  new_cmte_pot
end

# CLP is updated in place (i.e, cmte_indices updated), learn!() new cmte_pot and return that
function learn!(clp::SubsampleAppendCmteRetrain,
                  cmte_pot::CommitteePotential,
                  num_new_configs::Integer,
                  new_trainset::Vector{<:AbstractSystem})

  new_trainset_size = length(new_trainset)
  append_indices = (new_trainset_size-num_new_configs+1):new_trainset_size

  # append new indices (new configurations) to each cmte index set
  updated_indices = []
  for old_indices in clp.cmte_indices
    new_indices = reduce(vcat, [old_indices, [x for x in append_indices]])
    push!(updated_indices,new_indices)
  end

  clp.cmte_indices = updated_indices
  new_cmte_pot = learn(clp,cmte_pot,new_trainset)

  new_cmte_pot
end

# note that clp gets modified in place here
function retrain!(clp::SubsampleAppendCmteRetrain, sys::Molly.System, al::ALRoutine)
  new_trainset = al.trainset
  num_new_configs = length(al.cache[:trainset_changes]) # hard assumption that this is just a list of new systems
  @assert typeof(al.mlip) <: CommitteePotential
  old_cmte_pot = al.mlip

  new_cmte_pot = learn!(clp, old_cmte_pot, num_new_configs, new_trainset) # clp modified in place

  new_cmte_pot
end


#abstract type ALTriggerUpdate end

# This version can be used if I actually return updated updates
#function update_triggers(triggers, updates, sys::Molly.System, al::ALRoutine)
#  #new_trigger_updates = ALTriggerUpdate[]
#  new_trigger_updates = []
#  new_triggers = ActiveLearningTrigger[]
#
#  for (up,trigg) in zip(updates,triggers)
#    if !isnothing(up)
#      new_up, new_trigg = update_trigger(up, trigg;sys=sys,al=al)
#      push!(new_trigger_updates, new_up)
#      push!(new_triggers,new_trigg)
#    else
#      push!(new_trigger_updates, nothing)
#      push!(new_triggers,trigg)
#    end
#  end
#
#  return Tuple(new_trigger_updates), Tuple(new_triggers)
#end

#TODO is this the right order for the function arguments?
function update_triggers!(triggers, updates, sys::Molly.System, al::ALRoutine)
  new_triggers = ActiveLearningTrigger[]
  for (up,trigg) in zip(updates,triggers)
    if !isnothing(up)
      new_trigg = update_trigger!(up, trigg;sys=sys,al=al)
      push!(new_triggers,new_trigg)
    else
      push!(new_triggers,trigg)
    end
  end

  return Tuple(new_triggers)
end


function update_trigger!(update::SubsampleAppendCmteRetrain,
                        cmte_trigger::Union{CmteTrigger,SharedCmteTrigger};
                        al::ALRoutine,
                        kwargs...)

  old_cmte_pot = cmte_trigger.cmte_pot
  if isnothing(old_cmte_pot)
    @warn "Not actually updating CmteTrigger, assuming committee potential used for sys.general_inters and updated with retrain!()"
    return cmte_trigger
  else
    new_trainset = al.trainset
    num_new_configs = length(al.cache[:trainset_changes]) # hard assumption that this is just a list of new systems

    new_cmte_pot = learn!(update, old_cmte_pot, num_new_configs, new_trainset) # clp modified in place

    #is this kind of notation discouraged?
    new_cmte_trigger = typeof(cmte_trigger)(cmte_trigger; cmte_pot=new_cmte_pot)

    return new_cmte_trigger
  end
end

###===From makie/plot_contours.jl===####

## grid on 2d domain
function coord_grid_2d(
  limits::Vector{<:Vector},
  step::Real;
  dist_units = u"nm"
)
  xcoord = Vector(limits[1][1]:step:limits[1][2]) .* dist_units
  ycoord = Vector(limits[2][1]:step:limits[2][2]) .* dist_units
  return [xcoord, ycoord]
end

## generic potential energy function with coord as argument
function potential(inter, coord::SVector{2})
  sys = let coords=[coord]; () -> [SVector{2}(coords)]; end # pseudo-struct
  return AtomsCalculators.potential_energy(sys, inter)
end


## grid across potential energy surface below cutoff
function potential_grid_2d(
  inter,
  limits::Vector{<:Vector},
  step::Real;
  cutoff = nothing,
  dist_units = u"nm",
)
  rng1, rng2 = coord_grid_2d(limits, step; dist_units=dist_units)
  coords = SVector[]

  for i = 1:length(rng1)
      for j = 1:length(rng2)
          coord = SVector{2}([rng1[i],rng2[j]])
          Vij = ustrip(potential(inter, coord))
          if typeof(cutoff) <: Real && Vij <= cutoff
              append!(coords, [coord])
          elseif typeof(cutoff) <: Vector && cutoff[1] <= Vij <= cutoff[2]
              append!(coords, [coord])
          end
      end
  end

  return coords
end

###=========#####

abstract type IPErrors end

# train force and energy rmse, fisher divergence parameterized by integrator
struct Simple2DPotErrors <: IPErrors
  eval_int # if nothing, won't compute fisher divergence (e.g., for committee potential)
  compute_fisher::Bool
end

function initialize_error_metrics!(error_metric_type::Simple2DPotErrors, ddict::Dict)
  ddict["error_hist"] = Dict("rmse_e" => [],
                             "rmse_f" => [])
  if error_metric_type.compute_fisher
    ddict["error_hist"]["fd"] = []
  end
end

# like this is nearly identical to what she did but it's more verbose, so you have to explain why do it like this
function record_errors!(error_metric_type::Simple2DPotErrors, aldata::Dict, sys::Molly.System, al)
  r_e, r_f = compute_rmse(al.ref, al.mlip, error_metric_type.eval_int)
  append!(aldata["error_hist"]["rmse_e"], r_e)
  append!(aldata["error_hist"]["rmse_f"], r_f)
  if error_metric_type.compute_fisher
    fd = compute_fisher_div(al.ref, al.mlip, error_metric_type.eval_int)
    append!(aldata["error_hist"]["fd"], fd)
  end
end

struct DefaultALDataSpec
  error_metrics::IPErrors # Needs to be a stuct because of the initialization
  record_trigger_step::Bool
  #record_trigger_res::Bool
  record_parameters::Bool
  record_new_configs::Bool
  record_trigger_logs::Bool # could in theory only record a few triggers, but whatever
end

function DefaultALDataSpec(error_metrics::IPErrors;
          record_trigger_step::Bool=true,
          record_parameters::Bool=true,
          record_new_configs::Bool=false,
          record_trigger_logs::Bool=false)

  new_aldata_spec = DefaultALDataSpec(error_metrics,
                       record_trigger_step,
                       record_parameters,
                       record_new_config,
                       record_trigger_logs)
  new_aldata_spec
end

function initialize_al_record(al_spec::DefaultALDataSpec, al)
  aldata = Dict()
  initialize_error_metrics!(al_spec.error_metrics,aldata)
  if al_spec.record_trigger_step
    aldata["trigger_steps"] = []
  end

  if al_spec.record_parameters
    #this way the starting state is stored, but the "mlip_params" field will be the same length as the other fields
    aldata["original_mlip_params"] = get_params(al.mlip)
    aldata["mlip_params"] = []
  end

  if al_spec.record_new_configs
    aldata["new_configs"] = []
  end

  if al_spec.record_trigger_logs
    aldata["activated_trigger_logs"] = Dict()
    for trigger in al.triggers
      logger_ids = get_logger_ids(trigger) # can be more than one key
      filtered_logger_ids = filter(x->!isnothing(x), logger_ids)
      #setindex!.(Ref(aldata["activated_trigger_logs"]),
      #          [[] for _ in eachindex(filtered_logger_ids)],
      #          filtered_logger_ids) #this worked but ended up being more verbose than just the below for loop
      for logger_id in filtered_logger_ids
        aldata["activated_trigger_logs"][logger_id] = []
      end
    end
  end

  aldata
end


function record_al_record!(al_spec::DefaultALDataSpec, aldata::Dict, sys::Molly.System, al)
    record_errors!(al_spec.error_metrics, aldata, sys, al)

  if al_spec.record_trigger_step
    push!(aldata["trigger_steps"], al.cache[:trigger_step])
  end

  if al_spec.record_parameters
    #TODO need to formalize an interface to get parameters from mlips vs committee potentials
    params = get_params(al.mlip)
    push!(aldata["mlip_params"],params)
  end

  if al_spec.record_new_configs
    push!(aldata["new_configs"], al.cache[:trainset_changes])
  end

  if al_spec.record_trigger_logs
    for logger_key in keys(aldata["activated_trigger_logs"])
      observed_val = sys.loggers[logger_key].observable
      push!(aldata["activated_trigger_logs"][logger_key], observed_val)
    end
  end
end