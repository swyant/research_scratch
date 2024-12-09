import AtomsBase
using InteratomicPotentials
using ProgressBars
using LinearAlgebra: Diagonal, cond, norm

struct ChargedBasisPotential
    basis::InteratomicPotentials.BasisSystem
    α::Vector{Float64} #for pure basis 
    γ::Vector{Float64} #for basis*charge portion
    β::Vector{Float64} #for charge model
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

    β = AtA \ Atb

    β
end

function atomic_charges(A::AtomsBase.AbstractSystem, lbp::Union{ChargedBasisPotential,LBasisPotential}, Qtot::Float64=0.0; with_units=false)
    cld = compute_centered_descriptors(A,lbp.basis)    
    num_atoms = length(A)::Int64
    per_atom_Qtot = Qtot/num_atoms
    return atomic_charges(cld, lbp.β, per_atom_Qtot; with_units=with_units)
end

#function atomic_charges(cld, β, per_atom_Qtot::Float64=0.0; with_units=false)
#    cld = stack(cld)'
#
#    atomic_charges = cld*β
#
#    atomic_charges .+= per_atom_Qtot 
#
#    if with_units
#        return atomic_charges * u"e_au"
#    else
#        return atomic_charges
#    end
#end

function atomic_charges(cld, β, per_atom_Qtot::Float64=0.0; with_units=false)
    cld = stack(cld)'

    atomic_charges = cld*β

    #atomic_charges .+= per_atom_Qtot 
    full_atomic_charges = atomic_charges .+ per_atom_Qtot

    if with_units
        return full_atomic_charges * u"e_au"
    else
        return full_atomic_charges 
    end
end


function potential_energy1(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, Qtot::Float64=0.0)
    ld  = compute_local_descriptors(A,cbp.basis)
    cld = compute_centered_descriptors(A,cbp.basis)

    per_atom_Qtot = Qtot/length(A)

    qi  = atomic_charges(cld, cbp.β, per_atom_Qtot; with_units=false)
    
    # can't use eachindex
    pe=0.0
    for atm_idx in 1:length(A)
        pe += (cbp.α' * ld[atm_idx]) + cbp.γ'*(ld[atm_idx] * qi[atm_idx])
    end
    return pe
end 

function potential_energy2(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, Qtot::Float64=0.0)
    ld  = compute_local_descriptors(A,cbp.basis)
    cld = compute_centered_descriptors(A,cbp.basis)

    per_atom_Qtot = Qtot/length(A)

    qi  = atomic_charges(cld, cbp.β, per_atom_Qtot; with_units=false)
    
    pe = cbp.α' * sum(ld) + cbp.γ' * sum( ld .* qi)
    pe
end 

function potential_energy3(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, Qtot::Float64=0.0)
    ld  = compute_local_descriptors(A,cbp.basis)
    cld = compute_centered_descriptors(A,cbp.basis)

    per_atom_Qtot = Qtot/length(A)

    pe = cbp.α' * sum(ld)

    for i in 1:length(A)
        for k in 1:length(cbp.basis)
            for d in 1:length(cbp.basis)
                pe += cbp.γ[k]*cbp.β[d]*ld[i][k]*cld[i][d]
            end
         end
     end
    
    pe += per_atom_Qtot*cbp.γ'*sum(ld)    
    pe
end 

function potential_energy4(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, Qtot::Float64=0.0)
    
    basis_size = length(cbp.basis)

    ld  = compute_local_descriptors(A,cbp.basis)
    gd = sum(ld)
    cld = compute_centered_descriptors(A,cbp.basis)

    per_atom_Qtot = Qtot/length(A)

    pe = cbp.α' * gd 

    outer_ld_cld = zeros((basis_size, basis_size))

    for i in 1:length(A)
        outer_ld_cld += ld[i]*cld'[i]
    end

    pe += sum(cbp.γ*cbp.β' .* outer_ld_cld)
   
    pe += per_atom_Qtot*cbp.γ'*gd    
    pe
end 

# here the parameters are explicitly passed
function potential_energy5(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, ps::Vector{Float64}, Qtot::Float64=0.0)
    basis_size = length(cbp.basis)
    num_atoms = length(A)::Int64
    per_atom_Qtot = Qtot/length(A)

    α = ps[1:basis_size]
    β = ps[basis_size+1:2*basis_size]
    γ = ps[2*basis_size+1:end]

    ld = compute_local_descriptors(A,cbp.basis)
    gd = sum(ld)
    cld = stack(compute_centered_descriptors(A,cbp.basis))'

    outer_ld_cld = dropdims(sum(stack(ld[i]*cld[i,:]' for i in 1:num_atoms);dims=3),dims=3) #This should not be a one-liner

    pe = α'*gd + sum(γ*β' .* outer_ld_cld) + per_atom_Qtot*γ'*gd
    pe   
end


function potential_energy6(gd::Vector{Float64}, outer_ld_cld::Matrix{Float64}, cbp::ChargedBasisPotential, ps::Vector{Float64}, per_atom_Qtot::Float64=0.0)
    α = ps[1:basis_size]
    β = ps[basis_size+1:2*basis_size]
    γ = ps[2*basis_size+1:end]

    pe = α'*gd + sum(γ*β' .* outer_ld_cld) + per_atom_Qtot*γ'*gd
    pe   
end

function potential_energy7(gd::Vector{Float64}, 
                           outer_ld_cld::Matrix{Float64}, 
                           cbp::ChargedBasisPotential, 
                           α::Vector{Float64}, 
                           β::Vector{Float64}, 
                           γ::Vector{Float64},
                           per_atom_Qtot::Float64=0.0)

    pe = α'*gd + sum(γ*β' .* outer_ld_cld) + per_atom_Qtot*γ'*gd
    pe   
end

function potential_energy8(gd::Vector{Float64}, 
                           outer_ld_cld::Matrix{Float64}, 
                           cbp::ChargedBasisPotential, 
                           α, 
                           β, 
                           γ,
                           per_atom_Qtot::Float64=0.0)

    pe = α'*gd + sum(γ*β' .* outer_ld_cld) + per_atom_Qtot*γ'*gd
    pe   
end

function potential_energy9(gd::Vector{Float64}, 
                           outer_ld_cld::Matrix{Float64}, 
                           cbp::ChargedBasisPotential, 
                           β, 
                           γ,
                           per_atom_Qtot::Float64=0.0)

    pe = sum(γ*β' .* outer_ld_cld) + per_atom_Qtot*γ'*gd
    pe   
end

#ld is vector of vector 
#cld is matrix. 
#This is dumb
function compute_outer_ld_cld(ld, cld, num_atoms)
    basis_size = size(cld)[2]
    outer_ld_cld = zeros((basis_size,basis_size))
    for i in 1:num_atoms
        outer_ld_cld += ld[i] * cld[i,:]'
    end
    outer_ld_cld
end

function compute_outer_ld_cld2(ld, cld, num_atoms)
    outer_ld_cld = dropdims(sum(stack(ld[i]*cld[i,:]' for i in 1:num_atoms);dims=3),dims=3)
    outer_ld_cld
end

function prepare_train_db(configs::Vector{<:AtomsBase.FlexibleSystem}, basis::BasisSystem)
    dataset = Dict{String,Any}[]
    for config in configs
        num_atoms = length(config)
        ld = compute_local_descriptors(config,basis)
        gd = sum(ld)
        cld = compute_centered_descriptors(config,basis;ld=ld)
        _cld = stack(cld)'
        outer_ld_cld = compute_outer_ld_cld(ld,_cld,num_atoms)
        qtot = ustrip(get_total_charge(config))

        push!(dataset, Dict("config" => config,
                            "gd"     => gd,
                            "cld"    => cld,
                            "outer_ld_cld" => outer_ld_cld,
                            "qtot" => qtot,
                            "num_atoms" => num_atoms))
    end
    dataset        
end

function loss(ps, ds, cbp, qi_ref, eref; we=30.0, λ=0.01)
    basis_size = length(cbp.basis)
    α = ps[1:basis_size]
    β = ps[basis_size+1:2*basis_size]
    γ = ps[2*basis_size+1:end]

    loss = 0.0
    for (i,ddict) in enumerate(ds)
        qis_hat = atomic_charges(ddict["cld"],β,ddict["qtot"]/ddict["num_atoms"])
        loss += sum(.^(qis_hat .- qi_ref[i],2))
        
        e_hat = potential_energy8(ddict["gd"], 
                                  ddict["outer_ld_cld"],
                                  cbp, 
                                  α,
                                  β,
                                  γ,
                                  ddict["qtot"]/ddict["num_atoms"])

        loss += we/ddict["num_atoms"]*(e_hat-eref[i])^2
    end
    loss += λ*norm(ps)^2
    loss
end

function simpler_loss(ps, ds, cbp, qi_ref, eref; we=30.0, λ=0.01)
    basis_size = length(cbp.basis)
    β = ps[1:basis_size]
    γ = ps[basis_size+1:end]

    loss = 0.0
    for (i,ddict) in enumerate(ds)
        qis_hat = atomic_charges(ddict["cld"],β,ddict["qtot"]/ddict["num_atoms"])
        loss += sum(.^(qis_hat .- qi_ref[i],2))
        
        e_hat = potential_energy9(ddict["gd"], 
                                  ddict["outer_ld_cld"],
                                  cbp, 
                                  β,
                                  γ,
                                  ddict["qtot"]/ddict["num_atoms"])

        loss += we/ddict["num_atoms"]*(e_hat-eref[i])^2
    end
    loss += λ*norm(ps)^2
    loss
end

function ∇simpler_loss(ps, ds, cbp, qi_ref, eref; we=30.0, λ=0.01)
    basis_size = length(cbp.basis)
    β = ps[1:basis_size]
    γ = ps[basis_size+1:end]

    deriv = zero(ps)
    for (i,ddict) in enumerate(ds)
        outer_ld_cld = ddict["outer_ld_cld"]
        cld = stack(ddict["cld"])'
        qtot_per_atom = ddict["qtot"]/ddict["num_atoms"]
        gd = ddict["gd"]

        deriv[1:basis_size] += 2*cld'*cld*β - 2*cld'*qi_ref[i] # atomic charge loss contribution

        # derivative of potential energy
        pe_deriv = zero(ps)
        pe_deriv[1:basis_size] += outer_ld_cld'*γ
        pe_deriv[basis_size+1:end] += outer_ld_cld*β
        pe_deriv[basis_size+1:end] += qtot_per_atom*gd
        
        # energy loss derivative 
        e_hat = potential_energy9(gd,outer_ld_cld,cbp,β,γ,qtot_per_atom)
        pe_deriv *= 2*(we/ddict["num_atoms"])*(e_hat-eref[i])
        deriv += pe_deriv
    end

    deriv += 2*λ*ps
    deriv
end
#function check_derivative(sys_dict, cbp, ps)
#    basis_size = length(cbp.basis)
#    β = ps[1:basis_size]
#    γ = ps[basis_size+1:end]
#    e_hat = potential_energy9(sys_dict["gd"], 
#                              sys_dict["outer_ld_cld"],
#                              cbp, 
#                              β,
#                              γ,
#                              sys_dict["qtot"]/sys_dict["num_atoms"])
#    e_hat
#end

#function my_derivative(sys_dict,cbp,ps)
#    basis_size = length(cbp.basis)
#    β = ps[1:basis_size]
#    γ = ps[basis_size+1:end]
#
#    outer_ld_cld = sys_dict["outer_ld_cld"]
#    qtot_per_atom = sys_dict["qtot"]/sys_dict["num_atoms"]
#    gd = sys_dict["gd"]
#
#    deriv = zero(ps)
#    
#    deriv[1:basis_size] += outer_ld_cld'*γ
#    deriv[basis_size+1:end] += outer_ld_cld*β
#    deriv[basis_size+1:end] += qtot_per_atom*gd
#    
#    deriv
#end

#function check_derivative(sys_dict, cbp, ps,eref)
#    we = 30.0
#    basis_size = length(cbp.basis)
#    β = ps[1:basis_size]
#    γ = ps[basis_size+1:end]
#    e_hat = potential_energy9(sys_dict["gd"], 
#                              sys_dict["outer_ld_cld"],
#                              cbp, 
#                              β,
#                              γ,
#                              sys_dict["qtot"]/sys_dict["num_atoms"])
#    loss= we/sys_dict["num_atoms"]*(e_hat - eref)^2
#    loss
#end

# The ehat makes this slow, I suspect this could be optimized?
#function my_derivative(sys_dict,cbp,ps, eref)
#    basis_size = length(cbp.basis)
#    β = ps[1:basis_size]
#    γ = ps[basis_size+1:end]
#
#    outer_ld_cld = sys_dict["outer_ld_cld"]
#    qtot_per_atom = sys_dict["qtot"]/sys_dict["num_atoms"]
#    gd = sys_dict["gd"]
#
#    deriv = zero(ps)
#    
#    deriv[1:basis_size] += outer_ld_cld'*γ
#    deriv[basis_size+1:end] += outer_ld_cld*β
#    deriv[basis_size+1:end] += qtot_per_atom*gd
#    
#    e_hat = potential_energy9(gd,outer_ld_cld,cbp,β,γ,qtot_per_atom)
#    deriv = 2*(30.0/sys_dict["num_atoms"])*(e_hat-eref)*deriv
#
#    deriv   
#end


#function check_derivative(sys_dict, cbp, ps,qi_ref)
#    basis_size = length(cbp.basis)
#    β = ps[1:basis_size]
#
#    cld = sys_dict["cld"]
#    qtot_per_atom = sys_dict["qtot"]/sys_dict["num_atoms"]
#
#    qis_hat = atomic_charges(cld,β,qtot_per_atom)
#    loss = sum(.^(qis_hat-qi_ref,2))
#    loss 
#end
#
#function my_derivative(sys_dict,cbp,ps, qi_ref)
#    basis_size = length(cbp.basis)
#    β = ps[1:basis_size]
#
#    outer_ld_cld = sys_dict["outer_ld_cld"]
#    cld = stack(sys_dict["cld"])'
#    qtot_per_atom = sys_dict["qtot"]/sys_dict["num_atoms"]
#    gd = sys_dict["gd"]
#
#    deriv = zero(ps)
#    
#    deriv[1:basis_size] += 2*cld'*cld*β - 2*cld'*qi_ref
#    deriv   
#end

function check_derivative(configs, cbp, ps,eref, qi_ref)
    we = 30.0
    basis_size = length(cbp.basis)
    β = ps[1:basis_size]
    γ = ps[basis_size+1:end]
    loss = 0.0
    for (i,sys_dict) in enumerate(configs)
        qis_hat = atomic_charges(sys_dict["cld"],β,sys_dict["qtot"]/sys_dict["num_atoms"])
        loss += sum(.^(qis_hat .- qi_ref[i],2))

        e_hat = potential_energy9(sys_dict["gd"], 
                                  sys_dict["outer_ld_cld"],
                                  cbp, 
                                  β,
                                  γ,
                                  sys_dict["qtot"]/sys_dict["num_atoms"])
        loss += we/sys_dict["num_atoms"]*(e_hat - eref[i])^2
    end

    loss += 0.01*norm(ps)^2
    loss
end

# The ehat makes this slow, I suspect this could be optimized?
function my_derivative(configs,cbp,ps, eref, qi_ref)
    basis_size = length(cbp.basis)
    β = ps[1:basis_size]
    γ = ps[basis_size+1:end]
    deriv = zero(ps)
    
    for (i,sys_dict) in enumerate(configs)
        outer_ld_cld = sys_dict["outer_ld_cld"]
        qtot_per_atom = sys_dict["qtot"]/sys_dict["num_atoms"]
        cld = stack(sys_dict["cld"])'
        gd = sys_dict["gd"]

        deriv[1:basis_size] += 2*cld'*cld*β - 2*cld'*qi_ref[i]

        pe_deriv = zero(ps)
   
        pe_deriv[1:basis_size] += outer_ld_cld'*γ
        pe_deriv[basis_size+1:end] += outer_ld_cld*β
        pe_deriv[basis_size+1:end] += qtot_per_atom*gd
        e_hat = potential_energy9(gd,outer_ld_cld,cbp,β,γ,qtot_per_atom)
        pe_deriv *= 2*(30.0/sys_dict["num_atoms"])*(e_hat-eref[i])
        deriv += pe_deriv
    end
    deriv += 2*0.01*ps

    deriv   
end


function only_charge_loss(β, ds, qi_ref)
    loss = 0.0
    for (i,ddict) in enumerate(ds)
        qis_hat = atomic_charges(ddict["cld"],β,ddict["qtot"]/ddict["num_atoms"])
        loss += sum(.^(qis_hat .- qi_ref[i],2))
    end
    #loss += sqrt(0.01)*norm(β)
    loss += 0.01*norm(β)^2
    loss
end
