import AtomsBase
using InteratomicPotentials
using ProgressBars
using LinearAlgebra: Diagonal, cond

# -> Vector (size N) of Vectors (size K)
function compute_centered_descriptors(A::AtomsBase.AbstractSystem, basis::BasisSystem)
    num_atoms = length(A)::Int64

    ld = compute_local_descriptors(A,basis)
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

        if isapprox(0.0,total_charge;atol=1e-8)
            scale=1.0
        elseif isapprox(1.0,total_charge;atol=1e-8)
            scale=1.5
        elseif isapprox(-1.0, total_charge; atol=1e-8)
            scale=0.5
        else
            error("huh")
        end

        #A = scale*A
        A = exp(total_charge)*A
        
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

function atomic_charges(A::AtomsBase.AbstractSystem, lbp::LBasisPotential, Qtot::Float64=0.0; with_units=false)
    cld = compute_centered_descriptors(A,lbp.basis)    
    cld = stack(cld)'

    if isapprox(0.0,Qtot;atol=1e-8)
        scale=1.0
    elseif isapprox(1.0,Qtot;atol=1e-8)
        scale=1.5
    elseif isapprox(-1.0, Qtot; atol=1e-8)
        scale=0.5
    else
        error("huh")
    end

    #cld = scale*cld
    cld = exp(Qtot)*cld

    atomic_charges = cld*lbp.β

    num_atoms = length(A)::Int64
    atomic_charges .+= Qtot/num_atoms

    if with_units
        return atomic_charges * u"e_au"
    else
        return atomic_charges
    end
end
