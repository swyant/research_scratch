import AtomsBase
using InteratomicPotentials
using ProgressBars
using LinearAlgebra: Diagonal, cond

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

    β = AtA \ Atb

    β
end

function atomic_charges(A::AtomsBase.AbstractSystem, lbp::Union{ChargedBasisPotential,LBasisPotential}, Qtot::Float64=0.0; with_units=false)
    cld = compute_centered_descriptors(A,lbp.basis)    
    return atomic_charges(cld, lbp.β, qtot; with_units=with_units)
end

function atomic_charges(cld, β, Qtot::Float64=0.0; with_units=false)
    cld = stack(cld)'

    atomic_charges = cld*β

    num_atoms = length(A)::Int64
    atomic_charges .+= Qtot/num_atoms

    if with_units
        return atomic_charges * u"e_au"
    else
        return atomic_charges
    end
end

struct ChargedBasisPotential
    basis::InteratomicPotentials.BasisSystem
    γ::Vector{Float64} #for basis*charge portion
    β::Vector{Float64} #for charge model
end

function potential_energy(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, Qtot::Float64=0.0; ld=nothing, cld=nothing; qis=nothing)
    if isnothing(ld)
        ld = compute_local_descriptors(A,cbp.basis)
    end 

    if isnothing(cld)
        cld = compute_centered_descriptors(A, cbp.basis; ld=ld) 
    end 
    
    #if isnothing(qis)
    #    qis = atomic_charges(cld, cbp, Qtot) 
    #end

    pe = cbp.γ*(cld .* ld)
    pe
end

#This is probably not the right set of arguments, probably should dispatch off of cbp somehow...
function potential_energy(ld, cld, γ)
     pe = γ*(cld .* ld)
     pe 
end

#qi_ref <: Vector{Vector}
function loss(params,configs, qi_ref, eref)
    #compute residiuals of charges

    param_seglength = length(params)/2

    β = params[1:param_seglength]
    γ = paramas[param_seglength+1:end]
    
    loss = 0.0
    for (i,config) in enumerate(configs)
        qtot = get_total_charge(config)
        ld = compute_local_descriptors(config, cbp.basis)
        cld = compute_centered_descriptors(config, cbp.basis; ld=ld)

        qis_hat = atomic_charges(cld,β,qtot)
        loss += sum(.^(qis_hat .- qi_ref[i],2))
        
        e_hat = potential_energy(ld, cld, γ)
        loss += (e_hat - eref[i])^2
    end

    loss
end

res = optimize(loss, initial_guess)
