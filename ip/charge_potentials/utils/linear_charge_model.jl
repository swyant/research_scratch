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
    # should I pre-compute and store size of basis?
end

function Base.length(lcb::LinearChargeBasis)
    base_length = length(lcb.base_basis)
    #kappa_length = lcb.l_order*length(lcb.type_map)*base_length # with redundant zeros
    perelem_ld_length = base_length ÷ length(lcb.type_map) # to ensure it's integer
    kappa_length = lcb.l_order*length(lcb.type_map)*perelem_ld_length 
    total_length = base_length + kappa_length
end

function LinearChargeBasis(base_basis::BasisSystem,
                           charge_model::LBasisChargeModel,
                           l_order::Int64,
                           total_charge::Float64,
                           species_list::Vector{Symbol})
    type_map = Dict{Symbol, Int64}()
    for (i,species) in enumerate(species_list)
        type_map[species] = i
    end

    return LinearChargeBasis(base_basis, charge_model, l_order, total_charge, type_map)
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

    #electrostatic_descrs = compute_electrostatic_descriptors(A,lcb.type_map; qis=qis)

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

        #kappa terms (with redundant zeros)
        #kappa_term = zeros(Float64,num_elems*base_ld_size*lcb.l_order)
        #mixed_descriptor_mat = base_ld[i]' .* charge_descrs[i,:]
        #
        #kappa_start = (type_idx-1)*lcb.l_order*base_ld_size+1
        #kappa_end   = kappa_start + lcb.l_order*base_ld_size -1
        #kappa_term[kappa_start:kappa_end] .= reshape(mixed_descriptor_mat',:)        

        #new kappa terms
        perelem_ld_size = size(perelem_base_ld,2)
        kappa_term = zeros(Float64,num_elems*perelem_ld_size*lcb.l_order)
        mixed_descriptor_mat = perelem_base_ld[i,:]' .* charge_descrs[i,:]
        
        kappa_start = (type_idx-1)*lcb.l_order*perelem_ld_size+1
        kappa_end   = kappa_start + lcb.l_order*perelem_ld_size -1
        kappa_term[kappa_start:kappa_end] .= reshape(mixed_descriptor_mat',:)        


        #push!(ld, [alpha_term;gamma_term;kappa_term])
        push!(ld, [alpha_term;kappa_term])
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
function compute_electrostatic_descriptors(A,type_map,qis)
    num_elems = length(type_map)
    num_atoms = length(A)

    rijs = compute_rij(A)
    species = atomic_symbol(A)
    
    electro_descrs = [zeros(num_elems) for _ in 1:num_atoms]
    for i in 1:num_atoms
        qii = qis[i]
        for j in 1:num_atoms
            rij = rijs[i,j] # or r[j,i] whatever, symmetric
            tj = type_map[species[j]]
            electro_descrs[i][tj] += _fc(rij)*qii*qis[j]/rij
        end
    end
    electro_descrs
end

# Type stable?
#function compute_rij_mic(A)
#    cry_to_cart = Matrix(reduce(hcat, ustrip.(bounding_box(test_sys)) )')
#    cart_to_cry = inv(cry_to_cart)
#
#    cart_pos    =  reduce(hcat, ustrip.(position(A)))'
#    cry_pos     =  mapslices(x-> cart_to_cry*x, cart_pos, dims=2)
#
#    for 
#
#end

#function _compute_rij(a::Matrix{Float64})
#    n = size(a,1)
#    rij_mat = Matrix{Float64}(undef,n,n)
#    @inbounds for i in 1:n
#        rij_mat[i,i] =0.0
#        ai = view(a,i,:)
#        for j in (i+1):n
#            #rij_mat[i,j] = norm(_mic_adjust.(a[i,:] - a[j,:]))
#            rij_mat[i,j] = norm(a[i,:] - a[j,:])
#            rij_mat[j,i] = rij_mat[i,j]
#        end
#    end
#    rij_mat
#end

function _compute_rij(a::Matrix{Float64})
    a = permutedims(a,(2,1))
    n = size(a,2)
    rij_mat = Matrix{Float64}(undef,n,n)

    cry_delta_vec = zeros(3)
    cart_delta_vec = zeros(3)
     
    @inbounds for i in 1:n
        rij_mat[i,i] =0.0
        ai = view(a,:,i)
        for j in (i+1):n
            #rij_mat[i,j] = norm(_mic_adjust.(a[i,:] - a[j,:]))
            
            #rij_mat[i,j] = _my_eval(ai,view(a,:,j))
            
            #rij_mat[i,j] = norm(ai - view(a,:,j))
            
            #aj = view(a,:,j)
            #s = 0.0 
            #@inbounds for I in eachindex(ai,aj)
            #    aiI = ai[I]
            #    ajI = aj[I]
            #    s = s + abs2(aiI - ajI)
            #end
            #rij_mat[i,j] = sqrt(s)
            
            #aj = view(a,:,j)
            #rij_mat[i,j] = sqrt(sum(abs2.(ai .- aj)))
                      
            #aj = view(a,:,j)
            #rij_mat[i,j] = sqrt(sum(k -> abs2(ai[k] - aj[k]), eachindex(ai,aj)))           
            
            #aj = view(a,:,j)
            #tmp = [abs2(ai[k] - aj[k]) for k in eachindex(ai,aj)]
            #rij_mat[i,j] = sqrt(sum(tmp))           

            #aj = view(a,:,j)
            #for k in eachindex(ai,aj,cry_delta_vec)
            #    cry_delta_vec[k] = abs2(ai[k] - aj[k])
            #end
            #rij_mat[i,j] = sqrt(sum(cry_delta_vec))           
           
            aj = view(a,:,j)
            for k in eachindex(ai,aj,cry_delta_vec)
                cry_delta_vec[k] = abs2(_mic_adjust(ai[k] - aj[k]))
            end

            mul!(cart_delta_vec, cry_to_cart, cry_delta_vec) #unfortunately kills performance
            #cart_delta_vec = cry_to_cart*cry_delta_vec # But this is for sure worse
            rij_mat[i,j] = sqrt(sum(cart_delta_vec))           
            #rij_mat[i,j] = sqrt(sum(cry_delta_vec))           
          
            rij_mat[j,i] = rij_mat[i,j]
        end
    end
    rij_mat
end

function _my_eval(a,b)
    @inbounds begin 
    s = 0.0 #eval start
    for I in eachindex(a,b)
        aI = a[I]
        bI = b[I]
        s = _my_eval_reduce(s,_my_eval_op(aI,bI))
    end
    end
    return _my_eval_end(s)
end

function _my_eval_reduce(s,op_res)
    s + op_res
end

function _my_eval_op(a,b)
    s1 = abs(a-b)
    abs2(s1)
end

function _my_eval_end(s)
    sqrt(s)
end


function _mic_adjust(x)
    x > 0.5 ? x - 1.0 : (x < -0.5 ? x + 1.0 : x)
end

function _fc(rij)
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


#function compute_electrostatic_descriptors(A::AbstractSystem
