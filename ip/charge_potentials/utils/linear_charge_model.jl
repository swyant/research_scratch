using AtomsBase
import InteratomicPotentials: BasisSystem, compute_local_descriptors, potential_energy
using Unitful, UnitfulAtomic
using ProgressBars
using LinearAlgebra: Diagonal, cond, norm, mul!
# using SIMD 

#Probably should be parameterically typed
struct LBasisChargeModel
    basis::BasisSystem
    ξ::Vector{Float64}
end

#probably should be parametrically typed
# Assumes all systems have the same total charge
struct LinearChargeBasis <:BasisSystem
    base_basis::BasisSystem
    charge_model::LBasisChargeModel
    l_order::Int64
    total_charge::Float64 #TODO make all charge references integers
    type_map::Dict{Symbol,Int64}
    rcut::Float64
    # should I pre-compute and store size of basis?
end

function Base.length(lcb::LinearChargeBasis)
    base_length = length(lcb.base_basis)
    #kappa_length = lcb.l_order*length(lcb.type_map)*base_length # with redundant zeros
    perelem_ld_length = base_length ÷ length(lcb.type_map) # to ensure it's integer
    kappa_length = lcb.l_order*length(lcb.type_map)*perelem_ld_length 
    gamma_length = length(lcb.type_map)^2
    total_length = base_length + kappa_length + gamma_length
end

function LinearChargeBasis(base_basis::BasisSystem,
                           charge_model::LBasisChargeModel,
                           l_order::Int64,
                           total_charge::Float64,
                           species_list::Vector{Symbol},
                           rcut::Float64)
    type_map = Dict{Symbol, Int64}()
    for (i,species) in enumerate(species_list)
        type_map[species] = i
    end

    return LinearChargeBasis(base_basis, charge_model, l_order, total_charge, type_map, rcut)
end

# -> Vector (size N) of Vectors (size K)
function compute_centered_descriptors(A::AtomsBase.AbstractSystem, basis::BasisSystem; ld=nothing)
    
    num_atoms = length(A)::Int64
    
    if isnothing(ld)
        ld = compute_local_descriptors(A,basis)
    end
    mean_ld = sum(ld)/num_atoms
    centered_ld = ld .- (mean_ld,)

    centered_ld
end

# memory-efficient implementation
# No weights
function learn_charge_model(configs::Vector{<:AtomsBase.FlexibleSystem}, basis::BasisSystem; 
                            λ=0.01,
                            reg_style::Symbol=:default,
                            pbar=true)

    basis_size = length(basis)
    
    AtA = zeros(basis_size,basis_size)
    Atb = zeros(basis_size)

    if pbar
        iter = ProgressBar(configs)
    else
        iter = configs
    end

    for config in iter
        num_atoms = length(config)::Int64
        total_charge = ustrip(get_total_charge(config))

        b = atomic_charges = ustrip.(get_atomic_charges(config))
        b .-= total_charge/num_atoms

        A = centered_descrs = stack(compute_centered_descriptors(config,basis))'

        AtA .+= A'*A
        Atb .+= A'*b
    end 

    if reg_style == :default
        reg_matrix = λ*Diagonal(ones(size(AtA)[1]))
        AtA += reg_matrix
    elseif reg_style == :scale
        for i in 1:size(AtA,1)
           reg_elem = AtA[i,i]*(1+λ)
           reg_elem = max(λ,reg_elem)
           AtA[i,i] = reg_elem        
        end
    end

    println("condition number of AtA: $(cond(AtA))")

    ξ = AtA \ Atb

    ξ
end

function atomic_charges(A::AtomsBase.AbstractSystem, lbcm::LBasisChargeModel, Qtot::Float64=0.0; with_units=false)
    cld = compute_centered_descriptors(A,lbcm.basis)    
    num_atoms = length(A)::Int64
    per_atom_Qtot = Qtot/num_atoms
    return atomic_charges(cld, lbcm.ξ, per_atom_Qtot; with_units=with_units)
end


function atomic_charges(cld, ξ, per_atom_Qtot::Float64=0.0; with_units=false)
    cld = stack(cld)'

    atomic_charges = cld*ξ

    #atomic_charges .+= per_atom_Qtot 
    full_atomic_charges = atomic_charges .+ per_atom_Qtot # Zygote and mutable arrays

    if with_units
        return full_atomic_charges * u"e_au"
    else
        return full_atomic_charges 
    end
end

function compute_local_descriptors(A::AbstractSystem, lcb::LinearChargeBasis)
    Qtot = lcb.total_charge
    num_atoms = length(A)::Int64
    num_elems = length(lcb.type_map) 
    
    #TODO compute perelem ld instead of full ld vector, so no always-zero terms when computing kappa terms
    base_ld = compute_local_descriptors(A, lcb.base_basis) # Vector{Float64}
    base_ld_size = length(base_ld[1])

    #TODO if chargemodel and lcb don't share same basis, should not pass ld!
    @assert lcb.charge_model.basis == lcb.base_basis #probably a flag in the struct to indicate this?
    cld = compute_centered_descriptors(A, lcb.base_basis; ld=base_ld)
    
    # the initial one represents the beta term in Cuong's notation
    perelem_base_ld = hcat(ones(num_atoms),InteratomicPotentials.compute_perelem_local_descriptors(A,lcb.base_basis))

    qis = atomic_charges(cld, lcb.charge_model.ξ, Qtot/num_atoms)
    charge_descrs = compute_charge_descriptors(lcb.l_order, qis)

    electrostatic_descrs = compute_electrostatic_descriptors(A,lcb.type_map,qis,lcb.rcut)

    ld = Vector{Vector{Float64}}()
    species = atomic_symbol(A)
    for i in 1:num_atoms
        type_idx = lcb.type_map[species[i]]

        alpha_term = base_ld[i]

        ##gamma terms
        #gamma_term  = zeros(num_elems^2)
        #gamma_start = (type_idx-1)*num_elems+1
        #gamma_end   = gamma_start + num_elems -1 
        #gamma_term[gamma_start:gamma_end] .= electrostatic_descrs[i]

        gamma_term = electrostatic_descrs[i]

        #new kappa terms
        perelem_ld_size = size(perelem_base_ld,2)
        kappa_term = zeros(Float64,num_elems*perelem_ld_size*lcb.l_order)
        mixed_descriptor_mat = perelem_base_ld[i,:]' .* charge_descrs[i,:]
        
        kappa_start = (type_idx-1)*lcb.l_order*perelem_ld_size+1
        kappa_end   = kappa_start + lcb.l_order*perelem_ld_size -1
        kappa_term[kappa_start:kappa_end] .= reshape(mixed_descriptor_mat',:)        


        push!(ld, [alpha_term;gamma_term;kappa_term])
        #push!(ld, [alpha_term;kappa_term])
    end

    ld
end

function compute_charge_descriptors(l_order::Int64, qis::Vector{Float64})
    num_atoms = length(qis)
    charge_descrs = zeros((num_atoms,l_order))
    
    for l in 1:l_order
        charge_descrs[:,l] .= qis.^l
    end

    charge_descrs
end


# output array (length # of atoms) of arrays, where each inner array is length number of atoms
function compute_electrostatic_descriptors(A,type_map,qis,rcut)
    num_elems = length(type_map)
    num_atoms = length(A)

    rijs = compute_rij_mic(A)
    species = atomic_symbol(A)

    list_by_elem = [findall(==(sym), species) for sym in first.(sort(collect(type_map), by=last))]
    
    electro_descrs = [zeros(num_elems^2) for _ in 1:num_atoms]
    for i in 1:num_atoms
        qii = qis[i]
        ti = type_map[species[i]]
        # it's small so doesn't really matter, but allocating this array is not actually necessary
        perelem_electro_ld = zeros(num_elems) # nonzero components of full local descriptor for atom i 
        for eidx in 1:num_elems
            val = 0.0
            for j in list_by_elem[eidx]
                rij = rijs[i,j]
                if rij <= rcut && i != j# this is inefficient of course, would be better if this was already taken care of
                    val += qii*qis[j]*fc(rij,rcut)/rij
                end
            end
            perelem_electro_ld[eidx] = val
        end
        
        idx_start = (ti-1)*num_elems+1
        idx_end = idx_start + num_elems -1 

        electro_descrs[i][idx_start:idx_end] .= perelem_electro_ld
    end
    electro_descrs
end

# Type stable?, non-efficient due to cart -> crystal -> cart
# this won't be type stable, but the most expensive stuff is in theinner function which should be type stable
function compute_rij_mic(A::AtomsBase.AbstractSystem)
    cry_to_cart = Matrix(reduce(hcat, ustrip.(bounding_box(A)) )')
    cart_to_cry = inv(cry_to_cart)

    cart_pos    =  reduce(hcat, ustrip.(position(A)))'
    cry_pos     =  mapslices(x-> cart_to_cry*x, cart_pos, dims=2)

    rij_mat = _compute_rij_mic(cry_pos,cry_to_cart)
    rij_mat
end

function _compute_rij_mic(cry_pos, cry_to_cart::Matrix{Float64})
    cry_pos = permutedims(cry_pos,(2,1))
    n = size(cry_pos,2)
    rij_mat = Matrix{Float64}(undef,n,n)

    cry_delta_vec = zeros(3)
    cart_delta_vec = zeros(3)
     
    @inbounds for i in 1:n
        rij_mat[i,i] = 0.0
        ai = view(cry_pos,:,i)
        for j in (i+1):n           
            aj = view(cry_pos,:,j)
            for k in eachindex(ai,aj,cry_delta_vec)
                cry_delta_vec[k] = _mic_adjust(ai[k] - aj[k])
            end

            mul!(cart_delta_vec, cry_to_cart, cry_delta_vec) #unfortunately kills performance 
            rij_mat[i,j] = sqrt(sum(x-> abs2(x), cart_delta_vec))           
          
            rij_mat[j,i] = rij_mat[i,j]
        end
    end
    rij_mat
end

function _mic_adjust(x)
    x > 0.5 ? x - 1.0 : (x < -0.5 ? x + 1.0 : x)
end

function fc(rij, rcut)
    y = rij/rcut
    y2 = y^2

    y3 = 1.0 - y2*y;
    y4 = y3*y3 + 1e-6;
    y5 = sqrt(y4);
    y6 = exp(-1.0/y5);

    fcut = y6/exp(-1.0)
    fcut
end

function dfc(rij, rcut)
    y = rij/rcut
    y2 = y^2
    y3 = 1.0 - y2*y;
    y4 = y3*y3 + 1e-6;
    y5 = sqrt(y4);
    y6 = exp(-1.0/y5);
    y7 = y4*sqrt(y4)

    dfcut = ((3.0/(rcut*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7;
    return 1.0
end


# The only substantive difference between this and ooc_learn! is checking the charge (and it takes vector of configs, and the centered energies)
function learn_charge_energy_model(configs::Vector{<:AtomsBase.FlexibleSystem}, basis::BasisSystem, ref_energies::Vector{Float64};
                                λ=0.01,
                                reg_style::Symbol=:default,
                                pbar=true)

    basis_size = length(basis)
    
    AtA = zeros(basis_size,basis_size)
    Atb = zeros(basis_size)

    if pbar
        iter = ProgressBar(configs)
    else
        iter = configs
    end

    for (i,config) in enumerate(iter)
        ref_energy = ref_energies[i]

        #num_atoms = length(config)::Int64
        total_charge = ustrip(get_total_charge(config))
        if total_charge != basis.total_charge
            error("This system doesn't have the appropriate total charge!")
        end
        
        global_descrs = reshape(sum(compute_local_descriptors(config,basis)),:,1)'
        A = [global_descrs;]
        b = [ref_energy;]
        
        AtA .+= A'*A
        Atb .+= A'*b
    end 

    if reg_style == :default
        reg_matrix = λ*Diagonal(ones(size(AtA)[1]))
        AtA += reg_matrix
    elseif reg_style == :scale
        for i in 1:size(AtA,1)
           reg_elem = AtA[i,i]*(1+λ)
           reg_elem = max(λ,reg_elem)
           AtA[i,i] = reg_elem        
        end
    end

    println("condition number of AtA: $(cond(AtA))")

    ps = AtA \ Atb

    ps
end


############### Derivatives ##############  


#=
TODO How expensive is it to compute_force_descriptors and compute_peratom_force_descriptors?
could sum the latter to get the former, saves a call to POD library. Probably worth it
=#
function compute_charge_derivative_descriptors(config, lcb)
    num_atoms = length(config)::Int64
    global_force_descrs = compute_force_descriptors(config, lcb.basis)
    peratom_force_descrs = InteratomicPotentials.compute_peratom_force_descriptors(config, lcb.basis)

    # This is actually a bit faster than below, despice the repeated mapslices
    global_fd_mat = reduce(hcat, reduce(vcat, global_force_descrs)) # get rid of annoying vector{vector{vector}} format
    G = dropdims(mapslices(x-> lcb.ξ'*x,global_fd_mat, dims=1),dims=1)
    charge_derivs = dropdims(mapslices(x->lcb.ξ'*x,peratom_force_descrs,dims=3),dims=3)
    charge_derivs = mapslices(x-> x - G/num_atoms, charge_derivs, dims=1)

    #global_fd_mat = reduce(hcat, reduce(vcat, global_force_descrs))' # get rid of annoying vector{vector{vector}} format
    #centered_peratom_force_descrs = mapslices(x -> x- global_fd_mat ./ num_atoms, 
    #                                          peratom_force_descrs, 
    #                                          dims = (1,3))
    #charge_derivs = dropdims(mapslices(x->lcb.ξ'*x,centered_peratom_force_descrs,dims=3),dims=3)

    charge_derivs
end
