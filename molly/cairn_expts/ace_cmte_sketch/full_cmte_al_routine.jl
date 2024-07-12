import PotentialLearning: LinearProblem, SubsetSelector
import InteratomicPotentials: LinearBasisProblem
#=
One challenge I didn't execpt is that a simulator like NoseHoover basically has a different initialization routine (with more outputs) 
and in turn it's individual simulate_step!() has to track more

I might try to refactor Nose Hoover to have a cache field
=#

mutable struct ALRoutine
  ref 
  mlip
  trainset::Vector{<:AbstractSystem}
  triggers::Tuple{Varargs{<:ActiveLearningTrigger}}
  ss::SubsetSelector 
  lp::LearningProblem
  trigger_update::Union{Nothing,TriggerUpdate}
  sim_update::Union{Nothing,SimulatorUpdate}
  al_record_spec::ALRecordSpec
  cache::Dict
end


# initialize simulation
function initialize_sim!(
    sys::System,
    sim;
    n_threads::Integer=Threads.nthreads(),
    run_loggers=true,
    )

    sys.coords = wrap_coords.(sys.coords, (sys.boundary,))
    !iszero(sim.remove_CM_motion) && remove_CM_motion!(sys)
    nb = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    run_loggers!(sys, nb, 0, run_loggers; n_threads=n_threads)

    nb
end


function active_learn!(sys::Molly.System, 
                       sim::OverdampedLangevin, 
                       n_steps::Integer,
                       al::ALRoutine; 
                       n_threads::Integer=Threads.nthreads()
                       run_loggers=true,
                       rng=Random.GLOBAL_RNG)

  initialize_al_cache!(al)

  neighbors =  initialize_sim!(sys;, n_threads=n_threads, run_loggers=run_loggers)
  sys = initialize_triggers!(al.triggers,sys)
  aldata = initialize_al_record(al.al_record)

 
  for step_n in 1:n_steps
    neighbors = simulation_step!(sys,
                                 sim, 
                                 step_n,
                                 n_threads=n_threads, 
                                 neighbors=neighbors,
                                 run_loggers=run_loggers,
                                 rng=rng)

    # how to know which trigger activated? maybe return full result array
    #trigg_res, trigg_state = trigger_activated!(al.triggers,sys, step_n)
    trigg_res= trigger_activated!(al.triggers,sys, step_n)

    if trigg_res
      al.cache[:trigger_step] = step_n
      #al.cache[:trigger_state] = trigg_state
      
      # update trainset
      al.trainset, al.cache[:trainset_changes] = update_trainset!(al.ss,sys,al) # can modify caches
      # train new mlip
      al.mlip =  retrain!(al.lp,sys,al) # can modify caches
      sys.general_inters = (al.mlip,)

      if !isnothing(al.trigger_updates)
        trigger_updates, triggers = update_triggers(al, sys)
        al.trigger_updates = trigger_updates
        al.triggers = triggers
      end

      if !isnothing(al.sim_update)
        update_simulator!(al.sim_update,sim,sys,al)
      end

      record_al_record!(al.al_record,aldata,sys,al)
      reset_sys_cache_after_train(sys)
      reset_al_cache_after_train(al)
    end
    reset_sys_cache_after_step(sys)
  end
  aldata, sys
end

#=
Another tricky balance is how many arguments to pass. There's a case to be made of just offering all the information is the easiest, but then you risk not dispatching easily and type instabilities
=#

#struct WrapperLearningProblem <: LearningProblem
#  weights::Vector{<:Real}
#  intcpt::Boolean
#end 
#
## notice the extra boiler plate because I'm not just passing things in
#function retrain!(wlp::WrapperLearningProblem,sys::Molly.System,al)
#  current_lbp = al.mlip
#  trainset    = al.trainset
#
#  new_lbp = LinearBasisPotential(current_lbp.basis)
#  learn!(new_lbp,wlp.weights,wlp.intcpt)
#
#  new_lbp
#end

# Recomputing all descriptors, energies, forces for entire trainset
struct InefficientLearningProblem <: LearningProblem 
    weights::Vector{Float64}
    intcpt::Bool
end 

function InefficientLearningProblem(weights=[1000.0,1.0],intcpt=false)
    return InefficientLearningProblem(weights,intcpt)
end

function retrain!(ilp::InefficientLearningProblem, sys::Molly.System, al::AlRoutine) 
  lp = learn!(al.trainset, al.ref, al.mlip, ilp.weights, ilp.intcpt; e_flag=true, f_flag=true)
  new_mlip = deepcopy(al.mlip) #How to generalize? Should mlip be modified in place
  new_mlip.params = lp.β

  new_mlip
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

  new_trainset, new_sys
end

abstract struct 
  ss::SubsetSelector 
  lp::LearningProblem
  trigger_update::Union{Nothing,TriggerUpdate}
  sim_update::Union{Nothing,SimulatorUpdate}
  al_record_spec::ALRecordSpec
::AtomsCalcultors
