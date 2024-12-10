import AtomsBase
import InteratomicPotentials: potential_energy
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


function atomic_charges(cld, ζ, per_atom_Qtot::Float64=0.0; with_units=false)
    cld = stack(cld)'

    atomic_charges = cld*ζ

    #atomic_charges .+= per_atom_Qtot 
    full_atomic_charges = atomic_charges .+ per_atom_Qtot # Zygote and mutable arrays

    if with_units
        return full_atomic_charges * u"e_au"
    else
        return full_atomic_charges 
    end
end

function potential_energy_simple(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, ps, Qtot::Float64=0.0)
    ζ = ps[1:basis_size]
    α = ps[basis_size+1:2*basis_size]
    κ = ps[2*basis_size+1:end]
   
    ld  = compute_local_descriptors(A,cbp.basis)
    gd = sum(ld)
    cld = compute_centered_descriptors(A,cbp.basis)

    per_atom_Qtot = Qtot/length(A)

    qi  = atomic_charges(cld, ζ, per_atom_Qtot; with_units=false)
    
    pe = α'*gd + κ'*sum( ld .* qi)
    pe
end 


# here the parameters are explicitly passed
function potential_energy(A::AtomsBase.AbstractSystem, cbp::ChargedBasisPotential, ps::Vector{Float64}, Qtot::Float64=0.0)
    basis_size = length(cbp.basis)
    num_atoms = length(A)::Int64
    per_atom_Qtot = Qtot/length(A)

    ζ = ps[1:basis_size]
    α = ps[basis_size+1:2*basis_size]
    κ = ps[2*basis_size+1:end]

    ld = compute_local_descriptors(A,cbp.basis)
    gd = sum(ld)
    cld = stack(compute_centered_descriptors(A,cbp.basis))'

    outer_ld_cld = dropdims(sum(stack(ld[i]*cld[i,:]' for i in 1:num_atoms);dims=3),dims=3) #This should not be a one-liner

    pe = α'*gd + sum(κ*ζ' .* outer_ld_cld) + per_atom_Qtot*κ'*gd
    pe   
end


function potential_energy( gd::Vector{Float64}, 
                           outer_ld_cld::Matrix{Float64}, 
                           cbp::ChargedBasisPotential, 
                           ps,
                           per_atom_Qtot::Float64=0.0)
    basis_size = length(cbp.basis)
    ζ = ps[1:basis_size]
    α = ps[basis_size+1:2*basis_size]
    κ = ps[2*basis_size+1:end]
   
    pe = α'*gd + sum(κ*ζ' .* outer_ld_cld) + per_atom_Qtot*κ'*gd
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
    ζ = ps[1:basis_size]
    α = ps[basis_size+1:2*basis_size]
    κ = ps[2*basis_size+1:end]

    loss = 0.0
    for (i,ddict) in enumerate(ds)
        qis_hat = atomic_charges(ddict["cld"],ζ,ddict["qtot"]/ddict["num_atoms"])
        loss += sum(.^(qis_hat .- qi_ref[i],2))
        
        e_hat = potential_energy(ddict["gd"], 
                                 ddict["outer_ld_cld"],
                                 cbp, 
                                 ps,
                                 ddict["qtot"]/ddict["num_atoms"])

        loss += we/ddict["num_atoms"]*(e_hat-eref[i])^2
    end
    loss += λ*norm(ps)^2
    loss
end

function ∇loss(ps, ds, cbp, qi_ref, eref; we=30.0, λ=0.01)
    basis_size = length(cbp.basis)
    ζ = ps[1:basis_size]
    α = ps[basis_size+1:2*basis_size]
    κ = ps[2*basis_size+1:end]

    deriv = zero(ps)
    for (i,ddict) in enumerate(ds)
        outer_ld_cld = ddict["outer_ld_cld"]
        cld = stack(ddict["cld"])'
        qtot_per_atom = ddict["qtot"]/ddict["num_atoms"]
        gd = ddict["gd"]

        deriv[1:basis_size] += 2*cld'*cld*ζ - 2*cld'*qi_ref[i] # atomic charge loss contribution 2AtA - 2Ab

        # derivative of potential energy
        pe_deriv = zero(ps)
        pe_deriv[1:basis_size] += outer_ld_cld'*κ
        pe_deriv[basis_size+1:2*basis_size] += gd
        pe_deriv[2*basis_size+1:end] += outer_ld_cld*ζ
        pe_deriv[2*basis_size+1:end] += qtot_per_atom*gd
        
        # energy loss derivative 
        e_hat = potential_energy(gd,outer_ld_cld,cbp,ps,qtot_per_atom)
        pe_deriv *= 2*(we/ddict["num_atoms"])*(e_hat-eref[i])
        deriv += pe_deriv
    end

    deriv += 2*λ*ps
    deriv
end




function simpler_loss(ps, ds, cbp, qi_ref, eref; we=30.0, λ=0.01)
    basis_size = length(cbp.basis)
    ζ = ps[1:basis_size]
    κ = ps[basis_size+1:end]

    loss = 0.0
    for (i,ddict) in enumerate(ds)
        qis_hat = atomic_charges(ddict["cld"],ζ,ddict["qtot"]/ddict["num_atoms"])
        loss += sum(.^(qis_hat .- qi_ref[i],2))
        
        e_hat = potential_energy9(ddict["gd"], 
                                  ddict["outer_ld_cld"],
                                  cbp, 
                                  ζ,
                                  κ,
                                  ddict["qtot"]/ddict["num_atoms"])

        loss += we/ddict["num_atoms"]*(e_hat-eref[i])^2
    end
    loss += λ*norm(ps)^2
    loss
end

function ∇simpler_loss(ps, ds, cbp, qi_ref, eref; we=30.0, λ=0.01)
    basis_size = length(cbp.basis)
    ζ = ps[1:basis_size]
    κ = ps[basis_size+1:end]

    deriv = zero(ps)
    for (i,ddict) in enumerate(ds)
        outer_ld_cld = ddict["outer_ld_cld"]
        cld = stack(ddict["cld"])'
        qtot_per_atom = ddict["qtot"]/ddict["num_atoms"]
        gd = ddict["gd"]

        deriv[1:basis_size] += 2*cld'*cld*ζ - 2*cld'*qi_ref[i] # atomic charge loss contribution

        # derivative of potential energy
        pe_deriv = zero(ps)
        pe_deriv[1:basis_size] += outer_ld_cld'*κ
        pe_deriv[basis_size+1:end] += outer_ld_cld*ζ
        pe_deriv[basis_size+1:end] += qtot_per_atom*gd
        
        # energy loss derivative 
        e_hat = potential_energy9(gd,outer_ld_cld,cbp,ζ,κ,qtot_per_atom)
        pe_deriv *= 2*(we/ddict["num_atoms"])*(e_hat-eref[i])
        deriv += pe_deriv
    end

    deriv += 2*λ*ps
    deriv
end

function only_charge_loss(ζ, ds, qi_ref)
    loss = 0.0
    for (i,ddict) in enumerate(ds)
        qis_hat = atomic_charges(ddict["cld"],ζ,ddict["qtot"]/ddict["num_atoms"])
        loss += sum(.^(qis_hat .- qi_ref[i],2))
    end
    #loss += sqrt(0.01)*norm(ζ)
    loss += 0.01*norm(ζ)^2
    loss
end


function get_nonlinear_energies(ds,cbp,ps, all_qi_ref, all_e_ref)
    basis_size = length(cbp.basis)
    ζ = ps[1:basis_size]
    α = ps[basis_size+1:2*basis_size]
    κ = ps[2*basis_size+1:end]

    pred_energies = Float64[]
    pred_charges =  Vector{Float64}[]
    for config_dict in ds
        gd = config_dict["gd"]
        cld = config_dict["cld"]
        outer_ld_cld = config_dict["outer_ld_cld"]
        qtot_per_atom = config_dict["qtot"]/config_dict["num_atoms"]
        
        charges = atomic_charges(cld,ζ,qtot_per_atom)
        push!(pred_charges,charges)

        pe = potential_energy(gd,outer_ld_cld,cbp,ps,qtot_per_atom)
        push!(pred_energies,pe)
        
    end
    all_pred_charges = reduce(vcat, pred_charges)

    charge_mae,charge_rmse,charge_mape = calc_mae_rmse_mape(all_pred_charges,all_qi_ref)
    println("CHARGE\nmae: $(1000*charge_mae) me\nrmse: $(1000*charge_rmse) me\nmape: $(100*charge_mape) %\n")

    energy_mae,energy_rmse,energy_mape = calc_mae_rmse_mape(pred_energies,all_e_ref)
    println("ENERGY\nmae: $(1000*energy_mae) meV\nrmse: $(1000*energy_rmse) meV\nmape: $(100*energy_mape) %\n")
end
