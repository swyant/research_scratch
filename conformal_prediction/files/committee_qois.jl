export AbstractCommitteeQoI,
       CmteEnergy,
       CmteFlatForces,
       CmteForces,
       compute #probably too generic of a name to export

abstract type AbstractCommitteeQoI end

### Committee Energy QoIs ###

struct CmteEnergy <: AbstractCommitteeQoI
    cmte_reduce::Union{Nothing,Function}
    strip_units::Bool

    function CmteEnergy(reduce_function::Union{Nothing,Function}, strip_units::Bool)
        if !isnothing(reduce_function)
            if !_check_reduction_fn(reduce_function)
                error("Reduction function invalid. Must reduce to single Float, Integer, or Boolean")

            end
        end
        return new(reduce_function,strip_units)
    end

end

function CmteEnergy(reduce_fn::Union{Nothing, Function}=nothing;
                    strip_units::Bool=false)
    cmte_e_qoi = CmteEnergy(reduce_fn,strip_units)
    cmte_e_qoi
end


struct CmteEnergyCov <: AbstractCommitteeQoI
    strip_units::Bool
end

struct CmteAtomicEnergies <: AbstractCommitteeQoI
    cmte_reduce::Union{Nothing,Function}
    strip_units::Bool

    function CmteAtomicEnergies(reduce_function::Union{Nothing,Function}, strip_units::Bool)
        if !isnothing(reduce_function)
            if !_check_reduction_fn(reduce_function)
                error("Reduction function invalid. Must reduce to single Float, Integer, or Boolean")

            end
        end
        return new(reduce_function,strip_units)
    end
end

# Quickly hacking this up
function compute(qoi::CmteAtomicEnergies,config::PotentialLearning.Configuration,cmte_pot::CommitteePotential; cache_field=nothing)

    all_atomic_energies = compute_all_atomic_energies(config,cmte_pot)

    if qoi.strip_units
        if eltype(all_atomic_energies) <: Unitful.Quantity # may be unnecessary once AC interface firmly enforced?
            all_atomic_energies = ustrip.(all_atomic_energies)
        end
    end

    if !isnothing(qoi.cmte_reduce)
        all_atomic_energies = stack(all_atomic_energies) # num_atoms x num_members
        reduced_atomic_energies = qoi.cmte_reduce(all_atomic_energies,dims=2) # hard-coding which dimension to reduce over
        return reduced_atomic_energies
    else
        return all_atomic_energies
    end
end


# Quickly hacking this up
function compute(qoi::CmteAtomicEnergies,config::PotentialLearning.Configuration,cmte_pot::CommitteePotential; cache_field=nothing)

    all_atomic_energies = compute_all_atomic_energies(config,cmte_pot)

    if qoi.strip_units
        if eltype(all_energies) <: Unitful.Quantity # may be unnecessary once AC interface firmly enforced?
            all_energies = ustrip.(all_energies)
        end
    end

    if !isnothing(qoi.cmte_reduce)
        reduced_energies = qoi.cmte_reduce(all_energies)
        return reduced_energies
    else
        return all_energies
    end
end


function compute(qoi::CmteEnergy,sys::AbstractSystem,cmte_pot::CommitteePotential; cache_field=nothing)

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

    if qoi.strip_units
        if eltype(all_energies) <: Unitful.Quantity # may be unnecessary once AC interface firmly enforced?
            all_energies = ustrip.(all_energies)
        end
    end

    if !isnothing(qoi.cmte_reduce)
        reduced_energies = qoi.cmte_reduce(all_energies)
        return reduced_energies
    else
        return all_energies
    end
end

# Quickly hacking this up
function compute(qoi::CmteEnergy,config::PotentialLearning.Configuration,cmte_pot::CommitteePotential; cache_field=nothing)

    all_energies = compute_all_energies(config,cmte_pot)

    if qoi.strip_units
        if eltype(all_energies) <: Unitful.Quantity # may be unnecessary once AC interface firmly enforced?
            all_energies = ustrip.(all_energies)
        end
    end

    if !isnothing(qoi.cmte_reduce)
        reduced_energies = qoi.cmte_reduce(all_energies)
        return reduced_energies
    else
        return all_energies
    end
end

# Quickly hacking this up
function compute(qoi::CmteEnergyCov,config1::PotentialLearning.Configuration,config2::PotentialLearning.Configuration,
                cmte_pot::CommitteePotential; flip_second_sign=false, cache_field=nothing)

    all_energies1 = compute_all_energies(config1,cmte_pot)
    all_energies2 = compute_all_energies(config2,cmte_pot)

    if flip_second_sign 
        all_energies2 = -1 .* all_energies2
    end

    if qoi.strip_units
        if eltype(all_energies1) <: Unitful.Quantity # may be unnecessary once AC interface firmly enforced?
            all_energies1 = ustrip.(all_energies1)
            all_energies2 = ustrip.(all_energies2)
        end
    end

    return Statistics.cov(all_energies1, all_energies2)
end






### Committee Flattened Forces QoIs ###

struct CmteFlatForces <: AbstractCommitteeQoI
    cmte_reduce::Union{Nothing,Function}
    coord_and_atom_reduce::Union{Nothing,Function}
    reduce_order::Vector{Int64}
    strip_units::Bool

    function CmteFlatForces(cmte_reduce_function::Union{Nothing,Function},
                            coord_and_atom_reduce_function::Union{Nothing,Function},
                            reduce_order::Vector{Int64},
                            strip_units::Bool)
        num_reduce_fns = 0
        cmte_reduce_valid = true
        if !isnothing(cmte_reduce_function)
            num_reduce_fns += 1
            if !_check_reduction_fn(cmte_reduce_function)
                cmte_reduce_valid = false
            end
        end

        coord_and_atom_reduce_valid = true
        if !isnothing(coord_and_atom_reduce_function)
            num_reduce_fns += 1
            if !_check_reduction_fn(coord_and_atom_reduce_function)
                coord_and_atom_reduce_valid = false
            end
        end

        if !cmte_reduce_valid && coord_and_atom_reduce_valid
            error("Cmte reduction function invalid. Must reduce to single Float, Integer, or Boolean")
        elseif cmte_reduce_valid && !coord_and_atom_reduce_valid
            error("Coord_and_atom reduction function invalid. Must reduce to single Float, Integer, or Boolean")
        elseif !cmte_reduce_valid && !coord_and_atom_reduce_valid
            error("Cmte reduction and coord_and_atom reduction functions invalid. Must reduce to single Float, Integer, or Boolean")
        end

        if length(reduce_order) != num_reduce_fns
            error("Invalid reduce_order. Please use NamedTuple-based constructor")
        elseif length(reduce_order) > 0
            if maximum(reduce_order) > 2 || minimum(reduce_order) < 1
                error("Invalid reduce_order. Please use NamedTuple-based constructor")
            end
        end

        return new(cmte_reduce_function,coord_and_atom_reduce_function,reduce_order,strip_units)
    end

end

# There is some argument to be made that this NamedTuple-based constructor should be the inner constructor
function CmteFlatForces(nt::NamedTuple{<:Any, <:Tuple{Vararg{Function}}};
                        strip_units::Bool=false)
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
        error("""Only allowed keys are 'cmte', 'coord_and_atom' """)
    elseif length(nt) > 2
        error("There can be a maximum of 2 elements in the passed NamedTuple") # should never happen
    end

    reduce_order = Int64[]
    for (k,fn) in pairs(nt)
        push!(reduce_order,fn_dict[k]["idx"])
        fn_dict[k]["fn"] = fn
    end

    cmte_flat_force_qoi = CmteFlatForces(fn_dict[:cmte]["fn"],
                                         fn_dict[:coord_and_atom]["fn"],
                                         reduce_order,
                                         strip_units)
    cmte_flat_force_qoi
end

function CmteFlatForces(;strip_units=false)
    flat_force_qoi = CmteFlatForces(nothing,nothing,Int64[],strip_units)
    flat_force_qoi
end

function compute(qoi::CmteFlatForces,sys::AbstractSystem,cmte_pot::CommitteePotential; cache_field=nothing)
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

    if qoi.strip_units
        if eltype(all_flat_forces) <: Unitful.Quantity # may be unnecessary once AC interface firmly enforced?
            all_flat_forces = ustrip.(all_flat_forces)
        end
    end

    if isnothing(qoi.cmte_reduce) && isnothing(qoi.coord_and_atom_reduce)
        return all_flat_forces
    else
        inter = all_flat_forces
        for d in qoi.reduce_order
            inter = mapslices(reduce_dict[d],inter,dims=d)
        end

        if length(qoi.reduce_order) == 2
            @assert size(inter) == (1,1)
            final_qoi = inter[1]
        else
            for d in qoi.reduce_order
                @assert size(inter,d) == 1 # ensure singleton dimension
          end

          final_qoi = dropdims(inter,dims=Tuple(qoi.reduce_order))
        end
        return final_qoi
    end
end






### Committee Forces QoI ###

struct CmteForces <: AbstractCommitteeQoI
    cmte_reduce::Union{Nothing,Function}
    atom_reduce::Union{Nothing,Function}
    coord_reduce::Union{Nothing,Function}
    reduce_order::Vector{Int64}
    strip_units::Bool

    function CmteForces(cmte_reduce::Union{Nothing,Function},
                        atom_reduce::Union{Nothing,Function},
                        coord_reduce::Union{Nothing,Function},
                        reduce_order::Vector{Int64},
                        strip_units::Bool)

        # most of this extra logic is for more informative errors
        reduction_fns = [cmte_reduce,atom_reduce,coord_reduce]
        labels = ["cmte reduce", "atom reduce", "coord reduce"]
        reduction_invalid = [false,false,false]
        num_reduce_fns = 0

        for idx in eachindex(reduction_fns)
            if !isnothing(reduction_fns[idx])
                num_reduce_fns += 1
                if !_check_reduction_fn(reduction_fns[idx])
                    reduction_invalid[idx] = true
                end
            end
        end

        invalid_idxs = findall(reduction_invalid)
        if length(invalid_idxs) == 1
            e_msg = "The $(labels[invalid_idxs[1]]) function is invalid. Must reduce to single Float, Integer, or Boolean"
            error(e_msg)
        elseif length(invalid_idxs) > 1
            e_msg = "The "
            for iidx in invalid_idxs[1:end-1]
                e_msg = e_msg * "$(labels[iidx]), "
            end
            e_msg = e_msg * "$(labels[invalid_idxs[end]]) "
            e_msg = e_msg * "functions are invalid. Must reduce to single Float, Integer, or Boolean"
            error(e_msg)
        end

        if length(reduce_order) != num_reduce_fns
            error("Invalid reduce_order. Please use NamedTuple-based constructor")
        elseif length(reduce_order) > 0
            if maximum(reduce_order) > 3 || minimum(reduce_order) < 1
                error("Invalid reduce_order. Please use NamedTuple-based constructor")
            end
        end

        return new(cmte_reduce,atom_reduce,coord_reduce,reduce_order,strip_units)
    end


end

function CmteForces(;strip_units=false)
    cmte_force_qoi = CmteForces(nothing,nothing,nothing,Int64[],strip_units)
    cmte_force_qoi
end

function CmteForces(nt::NamedTuple{<:Any, <:Tuple{Vararg{Function}}}; strip_units=false)
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
        error("""Only allowed keys are 'cmte', 'atom', 'coord' """)
    elseif length(nt) > 3
        error("There can be a maximum of 3 elements in the passed NamedTuple") # should never happen
    end

    reduce_order = Int64[]
    for (k,fn) in pairs(nt)
        push!(reduce_order,fn_dict[k]["idx"])
        fn_dict[k]["fn"] = fn
    end

    cmte_force_qoi = CmteForces(fn_dict[:cmte]["fn"],
                                fn_dict[:atom]["fn"],
                                fn_dict[:coord]["fn"],
                                reduce_order,
                                strip_units)
    cmte_force_qoi
end

function compute(qoi::CmteForces,sys::AbstractSystem,cmte_pot::CommitteePotential; cache_field=nothing)

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

  if qoi.strip_units
      # Vector{Vector{<:AbstractArray{T,1}}} is the current format of raw_all_forces
      if eltype(raw_all_forces[1][1]) <: Unitful.Quantity # may be unnecessary once AC interface firmly enforced?
          raw_all_forces = [[ustrip.(forces_per_atom) for forces_per_atom in cmte_forces] for cmte_forces in raw_all_forces]
      end
  end

  if isnothing(qoi.cmte_reduce) && isnothing(qoi.atom_reduce) && isnothing(qoi.coord_reduce)
    return raw_all_forces
  else
    inter = stack(map(elem->stack(elem,dims=1),raw_all_forces), dims=1)  # num_cmte x num_atoms x 3
    for d in qoi.reduce_order
        inter = mapslices(reduce_dict[d],inter,dims=d)
    end

    if length(qoi.reduce_order) == 3
        @assert size(inter) == (1,1,1)
        final_qoi = inter[1]
    else
        for d in qoi.reduce_order
          @assert size(inter,d) == 1 # ensure singleton dimension
        end

        final_qoi = dropdims(inter,dims=Tuple(qoi.reduce_order))
    end
    return final_qoi
  end
end



### Utilities

# only run @ initialization of committee QoI
function _check_reduction_fn(fn::Function)
  test_arr1 = [1.0,2.0,3.0,4.0,5.0] # usually vector being operated on will be floats

  local res = nothing # to deal with the scoping in the try statement, can now use res in else
  try
    res = fn(test_arr1)
  catch e
    # note that a Unitful error is only really going to happen with the float array
    if isa(e,Unitful.DimensionError)
        expected_units = unit(e.x)
        unitted_test_arr = test_arr1 * expected_units
        try
            res = fn(unitted_test_arr)
        catch
        else
            if typeof(res) <: AbstractArray
                check_res = false
            else
                check_res = typeof(ustrip(res)) <: Union{<:Real, <:Integer, Bool}
            end
            return check_res
        end
    end
  else
    check_res = typeof(res) <: Union{<:Real, <:Integer, Bool}
    return check_res
  end

  test_arr2 = [1,2,3,4,5] # if prior reduction step produced integers
  try
    res = fn(test_arr2)
  catch
  else
    check_res = typeof(res) <: Union{<:Real, <:Integer, Bool}
    return check_res
  end

  test_arr3 = [true,true,false,true,false] # if prior reduction step produced integers
  try
    res = fn(test_arr3)
  catch
  else
    check_res = typeof(res) <: Union{<:Real, <:Integer, Bool}
    return check_res
  end

  false
end
