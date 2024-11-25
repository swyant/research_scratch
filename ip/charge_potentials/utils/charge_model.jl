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
        b = atomic_charges = ustrip.(get_atomic_charges(config))

        A = centered_descrs = stack(compute_centered_descriptors(config,basis))'

        AtA .+= A'*A
        Atb .+= A'*b
    end 

    reg_matrix = λ*Diagonal(ones(size(AtA)[1]))
    AtA += reg_matrix
    println("condition number of AtA: $(cond(AtA))")

    β = AtA \ Atb

    β
end

function atomic_charges(A::AtomsBase.AbstractSystem, lbp::LBasisPotential, Qtot=0.0; with_units=false)
    cld = compute_centered_descriptors(A,lbp.basis)    
    cld = stack(cld)'

    atomic_charges = cld*lbp.β

    if with_units
        return atomic_charges * u"e_au"
    else
        return atomic_charges
    end
end
