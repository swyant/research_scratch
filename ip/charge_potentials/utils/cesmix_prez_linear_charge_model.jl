using AtomsBase
import InteratomicPotentials: BasisSystem, compute_local_descriptors, potential_energy, compute_force_descriptors
using Unitful, UnitfulAtomic
using ProgressBars
using LinearAlgebra: Diagonal, cond, norm, mul!
# using SIMD

#Probably should be parameterically typed
# Also probably should include total charge here too? or only here actually
struct LBasisChargeModel
    basis::BasisSystem
    ξ::Vector{Float64}
end

struct LinearChargeBasis <:BasisSystem
    base_basis::BasisSystem
    charge_model::LBasisChargeModel
    l_order::Int64 #
    total_charge::Int64
    type_map::Dict{Symbol,Int64}
    rcut::Float64
    alpha::Bool
    kappa::Bool
    gamma::Bool
end

function Base.length(lcb::LinearChargeBasis)
    total_length = 0
    base_length = length(lcb.base_basis)
    if lcb.alpha
        total_length += base_length
    end
    #kappa_length = lcb.l_order*length(lcb.type_map)*base_length # with redundant zeros
    if lcb.gamma
        gamma_length = length(lcb.type_map)^2
        total_length += gamma_length
    end
    if lcb.kappa
        perelem_ld_length = base_length ÷ length(lcb.type_map) # to ensure it's integer
        kappa_length = lcb.l_order*length(lcb.type_map)*perelem_ld_length
        total_length += kappa_length
    end
    total_length
    #1
end


function LinearChargeBasis(base_basis::BasisSystem,
                           charge_model::LBasisChargeModel,
                           l_order::Int64,
                           total_charge::Float64,
                           species_list::Vector{Symbol},
                           rcut::Float64;
                           alpha=true,
                           kappa=true,
                           gamma=true)
    type_map = Dict{Symbol, Int64}()
    for (i,species) in enumerate(species_list)
        type_map[species] = i
    end

    return LinearChargeBasis(base_basis, charge_model, l_order, total_charge, type_map, rcut,alpha,kappa,gamma)
end

# -> Vector (size N) of Vectors (size K)
function compute_centered_descriptors(A::AtomsBase.AbstractSystem,
                                      basis::BasisSystem; ld=nothing)

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

    base_ld = compute_local_descriptors(A, lcb.base_basis) # Vector{Float64}
    base_ld_size = length(base_ld[1])
    cld = compute_centered_descriptors(A, lcb.base_basis; ld=base_ld)

    # the initial one represents the beta term in Cuong's notation
    perelem_base_ld = hcat(ones(num_atoms),
                          InteratomicPotentials.compute_perelem_local_descriptors(A,lcb.base_basis))
    if lcb.kappa || lcb.gamma
        qis = atomic_charges(cld, lcb.charge_model.ξ, Qtot/num_atoms)
        charge_descrs = compute_charge_descriptors(lcb.l_order, qis)
    end

    if lcb.gamma
        electrostatic_descrs = compute_electrostatic_descriptors(A,lcb.type_map,qis,lcb.rcut)
    end

    ld = Vector{Vector{Float64}}()
    species = atomic_symbol(A)
    for i in 1:num_atoms
        type_idx = lcb.type_map[species[i]]

        peratom_ld = Float64[]

        if lcb.alpha
            alpha_term = base_ld[i]
            peratom_ld = vcat(alpha_term,peratom_ld)
        end

        if lcb.gamma
            gamma_term = electrostatic_descrs[i]
            peratom_ld = vcat(peratom_ld,gamma_term)
        end

        if lcb.kappa
            #new kappa terms
            perelem_ld_size = size(perelem_base_ld,2)
            kappa_term = zeros(Float64,num_elems*perelem_ld_size*lcb.l_order)
            mixed_descriptor_mat = perelem_base_ld[i,:]' .* charge_descrs[i,:]

            kappa_start = (type_idx-1)*lcb.l_order*perelem_ld_size+1
            kappa_end   = kappa_start + lcb.l_order*perelem_ld_size -1
            kappa_term[kappa_start:kappa_end] .= reshape(mixed_descriptor_mat',:)
            peratom_ld = vcat(peratom_ld, kappa_term)
            #peratom_ld = vcat(peratom_ld, kappa_term[3])
        end

        push!(ld, peratom_ld)
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

    # this works, would probably be less bug prone if it was a dict
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

    r_mat, rij_mat = _compute_rij_mic(cry_pos,cry_to_cart)
    rij_mat
end

# Type stable?, non-efficient due to cart -> crystal -> cart
# this won't be type stable, but the most expensive stuff is in theinner function which should be type stable
function compute_rij_and_dispvec_mic(A::AtomsBase.AbstractSystem)
    cry_to_cart = Matrix(reduce(hcat, ustrip.(bounding_box(A)) )')
    cart_to_cry = inv(cry_to_cart)

    cart_pos    =  reduce(hcat, ustrip.(position(A)))'
    cry_pos     =  mapslices(x-> cart_to_cry*x, cart_pos, dims=2)

    r_mat, rij_mat = _compute_rij_mic(cry_pos,cry_to_cart)
    r_mat, rij_mat
end

function _compute_rij_mic(cry_pos, cry_to_cart::Matrix{Float64})
    cry_pos = permutedims(cry_pos,(2,1))
    n = size(cry_pos,2)
    rij_mat = Matrix{Float64}(undef,n,n)
    r_mat   = Array{Float64}(undef, n,n, 3)

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

            r_mat[i,j,:] = cart_delta_vec
            r_mat[j,i,:] = -1*cart_delta_vec
        end
    end
    r_mat, rij_mat
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
    return dfcut
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

# The only substantive difference between this and ooc_learn! is checking the charge (and it takes vector of configs, and the centered energies)
function learn_charge_energyforce_model(configs::Vector{<:AtomsBase.FlexibleSystem},
                                        basis::BasisSystem,
                                        ref_energies::Vector{Float64},
                                        ref_forces::Vector{Vector{Float64}};
                                        λ=0.01,
                                        ws = [30.0, 1.0],
                                        reg_style::Symbol=:default,
                                        pbar=true,
                                        eweight_normalized = :squared  # :squared, :standard or nothing
                                        )

    basis_size = length(basis)

    AtWA = zeros(basis_size,basis_size)
    AtWb = zeros(basis_size)

    if pbar
        iter = ProgressBar(configs)
    else
        iter = configs
    end

    for (i,config) in enumerate(iter)
        ref_energy = ref_energies[i]
        ref_force_comps = ref_forces[i]

        num_atoms = length(config)::Int64
        total_charge = ustrip(get_total_charge(config))
        if total_charge != basis.total_charge
            error("This system doesn't have the appropriate total charge!")
        end

        global_descrs = reshape(sum(compute_local_descriptors(config,basis)),:,1)'
        force_descrs  = stack(reduce(vcat,compute_force_descriptors(config,basis)))'

        A = [global_descrs; force_descrs]
        b = [ref_energy; ref_force_comps]


        if isnothing(eweight_normalized)
            we_norm = 1.0
        elseif eweight_normalized == :standard
            we_norm = 1/num_atoms
        elseif eweight_normalized == :squared
            we_norm = 1/num_atoms^2
        else
            error("eweight_normalized can only be nothing, :standard, or :squared")
        end

        W = Diagonal( [we_norm*ws[1];
                       ws[2]*ones(length(ref_force_comps))] )


        AtWA .+= A'*W*A
        AtWb .+= A'*W*b
    end

    if reg_style == :default
        reg_matrix = λ*Diagonal(ones(size(AtWA)[1]))
        AtWA += reg_matrix
    elseif reg_style == :scale
        for i in 1:size(AtWA,1)
           reg_elem = AtWA[i,i]*(1+λ)
           reg_elem = max(λ,reg_elem)
           AtWA[i,i] = reg_elem
        end
    end

    println("condition number of AtWA: $(cond(AtWA))")

    ps = AtWA \ AtWb

    ps
end


# The only substantive difference between this and ooc_learn! is checking the charge (and it takes vector of configs, and the centered energies)
function learn_energyforce_model(configs::Vector{<:AtomsBase.FlexibleSystem},
                                        basis::BasisSystem,
                                        ref_energies::Vector{Float64},
                                        ref_forces::Vector{Vector{Float64}};
                                        λ=0.01,
                                        ws = [30.0, 1.0],
                                        reg_style::Symbol=:default,
                                        pbar=true,
                                        eweight_normalized = :squared  # :squared, :standard or nothing
                                        )

    basis_size = length(basis)

    AtWA = zeros(basis_size,basis_size)
    AtWb = zeros(basis_size)

    if pbar
        iter = ProgressBar(configs)
    else
        iter = configs
    end

    for (i,config) in enumerate(iter)
        ref_energy = ref_energies[i]
        ref_force_comps = ref_forces[i]

        num_atoms = length(config)::Int64
        total_charge = ustrip(get_total_charge(config))

        global_descrs = reshape(sum(compute_local_descriptors(config,basis)),:,1)'
        force_descrs  = stack(reduce(vcat,compute_force_descriptors(config,basis)))'

        A = [global_descrs; force_descrs]
        b = [ref_energy; ref_force_comps]


        if isnothing(eweight_normalized)
            we_norm = 1.0
        elseif eweight_normalized == :standard
            we_norm = 1/num_atoms
        elseif eweight_normalized == :squared
            we_norm = 1/num_atoms^2
        else
            error("eweight_normalized can only be nothing, :standard, or :squared")
        end

        W = Diagonal( [we_norm*ws[1];
                       ws[2]*ones(length(ref_force_comps))] )


        AtWA .+= A'*W*A
        AtWb .+= A'*W*b
    end

    if reg_style == :default
        reg_matrix = λ*Diagonal(ones(size(AtWA)[1]))
        AtWA += reg_matrix
    elseif reg_style == :scale
        for i in 1:size(AtWA,1)
           reg_elem = AtWA[i,i]*(1+λ)
           reg_elem = max(λ,reg_elem)
           AtWA[i,i] = reg_elem
        end
    end

    println("condition number of AtWA: $(cond(AtWA))")

    ps = AtWA \ AtWb

    ps
end


############### Derivatives ##############


#=
TODO How expensive is it to compute_force_descriptors and compute_peratom_force_descriptors?
could sum the latter to get the former, saves a call to POD library. Probably worth it
=#
function compute_atomic_charge_derivatives(config, lcb)
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

    # descriptor derivatives (i.e., obtained from compute_force_descriptors, etc.)
    # have the negative built in on the assumption of it being for forces. So need to counteract that for charges
    -1.0 * charge_derivs
end

#= assuming that lcb.base_basis == lcb.charge_model.basis
OK so yes, this is very non-optimized.
The biggest saving will be the above described compute_peratom_force_descriptors --> global force descriptors --> used to compute base_global_force_descrs, atomic_charge_derivs, and peratom_fd.
=#
function compute_force_descriptors(A::AbstractSystem, lcb::LinearChargeBasis)
    num_atoms = length(A)::Int64
    num_elems = length(lcb.type_map)
    species = atomic_symbol(A)::Vector{Symbol}

    force_descrs =[ [Float64[] for alpha in 1:3] for _ in 1:num_atoms]

    if lcb.kappa || lcb.gamma
        base_ld = compute_local_descriptors(A, lcb.base_basis)
        @assert lcb.charge_model.basis == lcb.base_basis
        cld = compute_centered_descriptors(A, lcb.base_basis; ld=base_ld)
        qis = atomic_charges(cld, lcb.charge_model.ξ, lcb.total_charge/num_atoms)
        atomic_charge_derivs = compute_atomic_charge_derivatives(A,lcb.charge_model)
    end

    if lcb.kappa
        peratom_fd, peratom_fd_indices = InteratomicPotentials.compute_peratom_force_descriptors_withindices(A,lcb.base_basis)
        perelem_base_ld = hcat(ones(num_atoms), InteratomicPotentials.compute_perelem_local_descriptors(A, lcb.base_basis))
        perelem_ld_size = length(peratom_fd_indices[1])
        @assert perelem_ld_size == size(perelem_base_ld)[2]
        charge_descrs = compute_charge_descriptors(lcb.l_order, qis)

        reshaped_charge_descrs = reshape(charge_descrs,num_atoms,1,lcb.l_order)
        aug_charge_descrs = reshape(cat(ones(num_atoms), charge_descrs, dims=2)[:,1:end-1], num_atoms, 1, lcb.l_order)
        reshaped_perelem_base_ld = reshape(perelem_base_ld, num_atoms, perelem_ld_size,1)
        all_lorders = reshape([i for i in 1:lcb.l_order],1,1,lcb.l_order)
    end

    atom_types = [lcb.type_map[spec] for spec in species]
    peratom_fd_byelem = Dict{Int64, Array{Float64,3}}()
    elem_indices = Dict{Int64,Vector{Int64}}()
    num_atoms_byelem = Dict{Int64,Int64}()
    for eidx in 1:num_elems
        elem_atoms = findall(x-> x==eidx, atom_types)
        elem_indices[eidx] = elem_atoms
        num_atoms_byelem[eidx] = length(elem_atoms)
        if lcb.kappa && length(elem_atoms) > 0
            peratom_fd_byelem[eidx] = view(peratom_fd,:,elem_atoms,peratom_fd_indices[elem_atoms[1]]) # all atoms of the same type will have the same peratom_fd_indices, so just grab the first one
        end
    end


    #### Alpha term
    if lcb.alpha
        base_global_fd = compute_force_descriptors(A,lcb.base_basis)

        for j in 1:num_atoms
            for alpha in 1:3
                force_descrs[j][alpha] = vcat(force_descrs[j][alpha],base_global_fd[j][alpha])
            end
        end
    end


    #### Gamma term
    # This will need to be cleaned up eventually, need to precompute element indices by pairs, prefilter outer i == k, rik <=rcut
    # if you have the relevant indices, can do outer product (i.e., i, k) than sum over all the terms (i.e., sum over i and k)
    # then for every row_idx with atom j, can add in the fcut derivative selectively to relevant elem sets.
    if lcb.gamma
        disp_vecs, riks =  compute_rij_and_dispvec_mic(A)
        rcut = lcb.rcut
        for j in 1:num_atoms
            for alpha in 1:3
                row_idx = (j-1)*3 + alpha
                gamma_force_descrs = zeros(num_elems^2)
                for eidx1 in 1:num_elems
                    gamma_e1 = zeros(num_elems)
                    for i in elem_indices[eidx1]
                        dqi_dr = atomic_charge_derivs[row_idx, i]
                        qi = qis[i]
                        for eidx2 in 1:num_elems
                            for k in elem_indices[eidx2]
                                rik = riks[i,k]
                                if rik <= rcut && i != k
                                    dqk_dr = atomic_charge_derivs[row_idx,k]
                                    qk = qis[k]
                                    gamma_e1[eidx2] += (fc(rik,rcut)/rik)*(dqi_dr*qk + qi*dqk_dr)
                                    if j == i
                                        disp_vec_comp = disp_vecs[i,k,alpha]
                                        gamma_e1[eidx2] += qi*qk*disp_vec_comp * (rik*dfc(rik,rcut) - fc(rik,rcut))/rik^3
                                    elseif j == k
                                        disp_vec_comp = disp_vecs[i,k,alpha]
                                        gamma_e1[eidx2] -= qi*qk*disp_vec_comp * (rik*dfc(rik,rcut) - fc(rik,rcut))/rik^3
                                    end
                                end
                            end
                        end
                    end
                    idx_start = (eidx1-1)*num_elems + 1
                    idx_end = idx_start + num_elems -1
                    gamma_force_descrs[idx_start:idx_end] .= -1.0*gamma_e1 #no POD force descriptors here, have to manually account for -1
                 end
                force_descrs[j][alpha] = vcat(force_descrs[j][alpha], gamma_force_descrs)
            end
        end
    end


    ##### Kappa term
    if lcb.kappa
        @inbounds for j in 1:num_atoms
            for alpha in 1:3
                row_idx = (j-1)*3 + alpha
                kappa_force_descrs = zeros(perelem_ld_size*lcb.l_order*num_elems)
                for eidx in 1:num_elems
                    if num_atoms_byelem[eidx] == 0
                        continue
                    end
                    kappa_start = (eidx-1)*lcb.l_order*perelem_ld_size + 1
                    kappa_end = kappa_start + lcb.l_order*perelem_ld_size -1

                    term1 = reshape(view(peratom_fd_byelem[eidx],row_idx,:,:),num_atoms_byelem[eidx],perelem_ld_size,1) .* view(reshaped_charge_descrs,elem_indices[eidx],:,:)
                    term1 = dropdims(sum(term1,dims=1),dims=1)

                    term2 = view(reshaped_perelem_base_ld,elem_indices[eidx],:,:) .* all_lorders .* view(aug_charge_descrs,elem_indices[eidx],:,:) .* reshape(view(atomic_charge_derivs,row_idx,elem_indices[eidx]),num_atoms_byelem[eidx],1,1)
                    term2 = dropdims(sum(term2,dims=1),dims=1)

                    kappa_force_descrs[kappa_start:kappa_end] .= reshape(term1 .- term2, :)
                end
                force_descrs[j][alpha] = vcat(force_descrs[j][alpha], kappa_force_descrs)
             end
        end
    end

    force_descrs
end

