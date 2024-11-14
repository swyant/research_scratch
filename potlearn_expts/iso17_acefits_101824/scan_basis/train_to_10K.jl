using InteratomicPotentials, PotentialLearning
using Unitful
using Random
using JLD2
using LinearAlgebra: norm
using Statistics: mean 

include("../../md17_acefits_100124/splinter/pl_fit/subtract_peratom_e.jl")

function compute_error_metrics(ds_test, lb, vref_dict; 
                               e_test = get_all_energies(ds_test),
                               f_test = get_all_forces(ds_test))
    coeff_norm = norm(lb.β)
    @show coeff_norm

    e_test_pred = get_all_energies_w_onebody(ds_test, lb, vref_dict)
    e_test_mae, e_test_rmse, e_test_rsq = calc_metrics(e_test, e_test_pred)
    @show e_test_mae, e_test_rmse


    f_test_pred = get_all_forces(ds_test,lb)
    f_test_mae, f_test_rmse, f_test_rsq = calc_metrics(f_test_pred, f_test)
    @show f_test_mae, f_test_rmse
    
    return e_test_mae, e_test_rmse, f_test_mae, f_test_rmse, coeff_norm
end

ds_dict1 = load("10K_train_configs.jld2")
ds_dict2 = load("small_iso17_dataset.jld2")

base_ds_train = ds_dict1["ds_train"]
base_ds_test  = ds_dict2["ds_test"]

e_train = get_all_energies(base_ds_train)
e_test  = get_all_energies(base_ds_test)
f_test  = get_all_forces(base_ds_test)

num_atoms = length(get_system(base_ds_train[1])) # all systems have same number of atoms, entire ISO17 database

avg_energy_per_atom = mean(vcat(e_train,e_test))/num_atoms
vref_dict = Dict(:H => avg_energy_per_atom,
                 :C => avg_energy_per_atom,
                 :O => avg_energy_per_atom)

adjust_energies(base_ds_train,vref_dict)


# (body order, poly degree)
#ace_specs = [(3,12), (3,16), (3,20),
#             (4,10), (4,12), (4,14), (4,15), (4,16),
#             (5,9), (5,10), (5,11), (5,12)]

ace_specs = [(4,15), (4,16), (5,12)]


ddict = Dict()
for spec in ace_specs
    println("body order : $(spec[1])  poly degree: $(spec[2])")
    ace = ACE(species           = [:C,:O,:H],
              body_order        = spec[1],
              polynomial_degree = spec[2],
              wL                = 2.0,
              csp               = 1.0,
              r0                = 1.43,
              rcutoff           = 4.4 )
    println("length of ace: $(length(ace))") 
      
    lb = LBasisPotential(ace)
    ###### fitting procedure
    println("ooc fitting")
    ws = [30.0, 1.0]
    _AtWA, _AtWb = PotentialLearning.ooc_learn!(lb, base_ds_train;ws=ws,symmetrize=false, λ=0.01)

    println("computing test descriptors")
    edescr_test = compute_local_descriptors(base_ds_test, ace)
    fdescr_test = compute_force_descriptors(base_ds_test, ace)
    ds_test = DataSet(base_ds_test .+ edescr_test .+ fdescr_test)


    e_test_mae, e_test_rmse, f_test_mae, f_test_rmse, coeff_norm  = compute_error_metrics(
                                                                        ds_test,
                                                                        lb,
                                                                        vref_dict;
                                                                        e_test = e_test,
                                                                        f_test = f_test)
    ddict[spec] = Dict("coeffs" => lb.β,
                       "length" => length(ace),                        
                       "norm"   => coeff_norm,
                       "e_mae"  => e_test_mae,
                       "e_rmse" => e_test_rmse,
                       "f_mae"  => f_test_mae,
                       "f_rmse" => f_test_rmse)

    # I just need to be very cautious with memory usage                                                                                             
    ds_test      = nothing
    edescr_test  = nothing
    fdescr_test  = nothing  
    lb           = nothing 
    _AtWA        = nothing
    _AtWb        = nothing

    GC.gc()
end

save("summary_of_ace_basis_scan_iso17_fit-to-10K.jld2", Dict("ddict" => ddict))
