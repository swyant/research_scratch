import PotentialLearning: get_random_subset, get_inclusion_prob, KernelMatrix


function KernelMatrix(
    F::Union{Vector{Vector{T}},Vector{Symmetric{T,Matrix{T}}}},
    k::Kernel,
) where {T}
    n = length(F)
    K = zeros(n, n)
    for i = 1:n
        for j = i:n
            K[i, j] = PotentialLearning.compute_kernel(F[i], F[j], k)
        end
    end
    Symmetric(K)
end

## compute ell ensemble for kDPP
function compute_ell_ensemble(
    ds::DataSet,
    f::Feature,
    k::Kernel;
    dt = LocalDescriptors,
)
    K = KernelMatrix(ds, f, k; dt = dt)
    ell = EllEnsemble(K)
    return K, ell
end

function compute_ell_ensemble(
    descriptors::Union{Vector{Vector{T}},Vector{Symmetric{T,Matrix{T}}}},
    k::Kernel
) where T <: Real
    K = KernelMatrix(descriptors, k)
    ell = EllEnsemble(K)
    return K, ell
end


## get random subset by kDPP
function get_random_subset(
    ell::EllEnsemble{Float64},
    batch_size::Int64,
)
    # compute kDPP
    rescale!(ell, batch_size)
    dpp = kDPP(ell, batch_size)
    # draw subset
    inds = get_random_subset(dpp)   

    return inds
end


## 
function get_inclusion_prob(
    ell::EllEnsemble{Float64},
    batch_size::Int64
)
    # compute kDPP
    rescale!(ell, batch_size)
    dpp = kDPP(ell, batch_size)
    # compute inclusion probabilities
    return get_inclusion_prob(dpp)
end


