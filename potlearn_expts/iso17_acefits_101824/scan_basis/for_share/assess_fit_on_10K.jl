using InteratomicPotentials, PotentialLearning
using Unitful
using Random
using JLD2, CSV, DataFrames
using LinearAlgebra: norm
using Statistics: mean 

include("../../../md17_acefits_100124/splinter/pl_fit/subtract_peratom_e.jl")

test_xyzfile  = "../../../../../../datasets/ISO17/my_iso17_test.xyz"
raw_test   = load_data(test_xyzfile, ExtXYZ(u"eV", u"Å"))
test_idxs  = vec(Matrix(CSV.read("../../initial_exploration/10K_random_test_idxs",DataFrame,header=false)))
ds_test    = raw_test[test_idxs]


ds1_dict = load("../small_iso17_dataset.jld2")
ds_test_for_vref  = ds1_dict["ds_test"]

ds2_dict = load("../10K_train_configs.jld2")
ds_train = ds2_dict["ds_train"]

# computing average energy, to subtract from the dataset
e_train = get_all_energies(ds_train)
e_test_for_vref  = get_all_energies(ds_test_for_vref)
e_Test  = get_all_energies(ds_test)
f_test  = get_all_forces(ds_test) #just computing this here, used later for error metrics

num_atoms = length(get_system(ds_train[1])) # all systems have same number of atoms for entire ISO17 database
avg_energy_per_atom = mean(vcat(e_train,e_test_for_vref))/num_atoms
vref_dict = Dict(:H => avg_energy_per_atom,
                 :C => avg_energy_per_atom,
                 :O => avg_energy_per_atom)


ace = ACE(species           = [:C,:O,:H],
          body_order        = 4,
          polynomial_degree = 16,
          wL                = 2.0,
          csp               = 1.0,
          r0                = 1.43,
          rcutoff           = 4.4 )
println("length of ace: $(length(ace))") 
  
lb = LBasisPotential(ace)

fitted_params = vec(Matrix(CSV.read("./4body_16poly_fit-to-10K_beta.txt",DataFrame,header=false)))

# Computing the test error
println("computing test descriptors")
@time begin
edescr_test = compute_local_descriptors(ds_test, ace)
fdescr_test = compute_force_descriptors(ds_test, ace)
ds_test = DataSet(ds_test .+ edescr_test .+ fdescr_test)
end

coeff_norm = norm(lb.β)
@show coeff_norm

# modified get all_energies, basiically adding back the one-body energy terms
println("running test predictions")
@time begin
e_test_pred = get_all_energies_w_onebody(ds_test, lb, vref_dict) 
e_test_mae, e_test_rmse, e_test_rsq = calc_metrics(e_test, e_test_pred)
@show e_test_mae

f_test_pred = get_all_forces(ds_test,lb)
f_test_mae, f_test_rmse, f_test_rsq = calc_metrics(f_test_pred, f_test)
@show f_test_mae
end

save("10K_test_descrs.jld2", Dict("edescr_test" => edescr_test, 
                                  "fdescr_test" => fdescr_test))
