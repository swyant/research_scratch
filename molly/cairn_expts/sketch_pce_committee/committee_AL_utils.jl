using Random: randperm
using Statistics: std
using  Molly: Molly, System

# As it is, this struct has to carry around the fit indices
# which is not broadly applicable to all committee_potentials
# Will need a better soln eventually (probably lumped in with a training struct to pass to AlRoutine)
struct Committee_Potential 
  pot_members::Vector{MLInteraction} #Should be AbstractPotential
  leader::Integer
  fit_idxs::Vector{Vector{Int64}} 
end 

function Committee_Potential(members::Vector{MLInteraction}, fit_idxs::Vector{Vector{Int64}})
  cmte_pot = Committee_Potential(members,1,fit_idxs)
  cmte_pot
 end 

#ens_stat_func 
struct Committee_Trigger <: ActiveLearningTrigger
    ens_stat_func::Function 
    eval::Function 
    thresh::Real
end


# This will allow for an "external" committe_potential
function MaxCmteEStdev(; thresh::Real=0.1, cmte_pot::Committee_Potential=nothing)
  if !isnothing(cmte_pot)
      _estdev(; kwargs...) = _energy_stdev(; cmte_pot=cmte_pot, kwargs...)
      return LowerThreshold(_estdev, thresh)
  else
    return LowerThreshold(energy_stdev, thresh)
end 


function _energy_stdev(; sys::Molly.System, cmte_pot::Committee_Potential=nothing)
  if isnothing(cmte_pot) # assume the general_inters of the system is a committee potential
      cmte_pot = sys.general_inters[1]
  stdev_e = energy_stdev(sys,cmte_pot)
  stdev_e
end


# The "natural" function that takes in a system and committee potential, outputs energy stdev
# maybe could generalize to any summary statistic?
function energy_stdev(sys::AbstractSystem, cmte_pot::Committee_Potential=nothing) # reverse the arguments?
  energies = []
  #energies = all_potential_energies(cmte_pot,sys) # could dispatch here instead, probably better
  for pot in cmte_pot.pot_members
     potential_energy(pce,sys) 
  end

  stdev_e = std(energies)
  stdev_e
end



function energies_statistic(sys, cmte_pot, statistic=default)

#=
eval(; kwargs...) < thresh ==> energy_stdev(; sys=sys, cmte_pot=cmte_pot)
eval(; kwargs...) < thresh ==> _energy_stdev(; sys=sys, cmte_pot=cmte_pot) ==> energy_stdev(sys,cmte_pot) *my preferred soln
eval(; kwargs...) < thresh ==> energy_stdev(; sys=sys, cmte_pot=cmte_pot) ==> _energy_stdev(sys,cmte_pot)
eval(sys, uq_metric; kwargs...) < thresh ==>  (sys,uq_metric) -> energy_stdev(sys,uq_metric.cmte_pot)


I do wonder, if part of what I'm trying to do is to ensure that trigger_activated() in active_learn!() has a fixed calling pattern.
But maybe we jsut decide that trigger_activated takes in the trigger and the system?
Otherwise we can have it receive a named tuple, where that is previously generated in lines above based off the trigger type and (other fields?)
=#
#= 
Why do it this way: I may want to naively take a random subset of the full data, or 
 I already have a list of subset ind, and I'm taking an additional subset of that, but 
 I still want the ind to be wrt the original dataset
=#
function obtain_train_inds(frac; curr_subset_inds=nothing, set_size=nothing)
    @assert frac <= 1.0
    if !isnothing(curr_subset_inds)
        if !isnothing(set_size)
            print("set size will be overwritten")
        end
        set_size = length(curr_subset_inds)
    end
    num_select = Int(floor(frac*set_size))

    perm_inds = randperm(set_size)
    rand_set_inds = perm_inds[begin:1:num_select]

    if isnothing(curr_subset_inds)
        return rand_set_inds
    else 
        return curr_subset_inds[rand_set_inds]
    end
end


function train_committee!(template_pot::MLInteraction, ref_pot, full_train_set::Vector{T}, num_members::Integer; frac=0.5) where T<:AbstractSystem
  pot_members = Vector{MLInteraction}()
  all_train_inds = Vector{Vector{Int64}}()
  for i in 1:num_members
      train_inds = obtain_train_inds(frac; set_size=length(full_train_set))
      push!(all_train_inds,train_inds)

      pot = deepcopy(template_pot)
      pot = train_potential_e!(full_train_set[train_inds],ref_pot,pot)
      push!(pot_members,pot)
  end

  cmte_pot = committee_potential(pot_members,all_train_inds)

  cmte_pot, all_train_inds  
end


function my_active_learn!(sys::System,
    sim::OverdampedLangevin,
    n_steps::Integer,
    al::ActiveLearnRoutine;
    n_threads::Integer=Threads.nthreads(),
    run_loggers=true,
    rng=Random.GLOBAL_RNG,
)
    sys.coords = wrap_coords.(sys.coords, (sys.boundary,))
    !iszero(sim.remove_CM_motion) && remove_CM_motion!(sys)
    neighbors = find_neighbors(sys, sys.neighbor_finder; n_threads=n_threads)
    run_loggers!(sys, neighbors, 0, run_loggers; n_threads=n_threads)
    compute_error_metrics!(al)

    ct = 0
    for step_n in 1:n_steps
        ct += 1
        neighbors = simulation_step!(sys,
                        sim,
                        step_n,
                        n_threads=n_threads,
                        neighbors=neighbors,
                        run_loggers=run_loggers,
                        rng=rng,
        )

        # online active learning
        if trigger_activated(al.trigger; ens_old=al.sys_train, sys_new=sys, step_n=step_n) && ct >= al.burnin
            println("train on step $step_n")
            al.sys_train = al.update_func(sim, sys, al.sys_train)
            al.train_func(sys, al.sys_train, al.ref) # retrain potential
            al.mlip = sys.general_inters[1]
            append!(al.train_steps, step_n)
            append!(al.param_hist, [sys.general_inters[1].params])
            compute_error_metrics!(al)
            ct = 0 # reset counter
        end
    end
    return al
end


#= 
Issues with the approach: 

I have to use the committee potential as the driver of the 
=#
