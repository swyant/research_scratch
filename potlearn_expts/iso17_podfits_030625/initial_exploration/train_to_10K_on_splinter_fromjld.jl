using InteratomicPotentials, PotentialLearning
using Unitful
using Random
using JLD2, CSV, DataFrames
using LinearAlgebra: norm
using Statistics: mean 

include("../files/subtract_peratom_e.jl")

#train_xyzfile = "../../../../../../datasets/ISO17/my_iso17_train.xyz"
#test_xyzfile  = "../../../../../../datasets/ISO17/my_iso17_test.xyz"
#
##train_idxs = Random.randperm(length(raw_train))[1:10000]
##test_idxs  = Random.randperm(length(raw_test))[1:10000]
##CSV.write("10K_random_train_idxs", DataFrame(Tables.table(train_idxs)),header=false)
##CSV.write("10K_random_test_idxs", DataFrame(Tables.table(test_idxs)),header=false)
#println("loading training xyz")
#@time begin
#raw_train  = load_data(train_xyzfile, ExtXYZ(u"eV", u"Å"))
#train_idxs = vec(Matrix(CSV.read("../10K_random_train_idxs",DataFrame,header=false)))
#ds_train   = raw_train[train_idxs]
#end
#
#println("loading test xyz")
#@time begin
## only taking 500 test configs. Note these are *unknown/unseen* molecules
#raw_test   = load_data(test_xyzfile, ExtXYZ(u"eV", u"Å"))
#test_idxs  = vec(Matrix(CSV.read("../../initial_exploration/10K_random_test_idxs",DataFrame,header=false)))[1:500]
#ds_test    = raw_test[test_idxs]
#end
#

ds1_dict = load("../files/small_iso17_dataset.jld2")
ds_test  = ds1_dict["ds_test"]

ds2_dict = load("../files/10K_train_configs.jld2")
ds_train = ds2_dict["ds_train"]

# computing average energy, to subtract from the dataset
e_train = get_all_energies(ds_train)
e_test  = get_all_energies(ds_test)
f_test  = get_all_forces(ds_test) #just computing this here, used later for error metrics

num_atoms = length(get_system(ds_train[1])) # all systems have same number of atoms for entire ISO17 database
avg_energy_per_atom = mean(vcat(e_train,e_test))/num_atoms
vref_dict = Dict(:H => avg_energy_per_atom,
                 :C => avg_energy_per_atom,
                 :O => avg_energy_per_atom)

# only do this for the train set, not the test set
adjust_energies(ds_train,vref_dict) # This permanently changes the energies in the dataset

param_file = "4body_COH_param.pod"
lmp_pod = LAMMPS_POD(param_file, [:C, :O, :H])

println("length of POD basis: $(length(lmp_pod))") 
  
lb = LBasisPotential(lmp_pod)
###### fitting procedure
println("ooc fitting")
ws = [30.0, 1.0]

# no intercept, Train force/energy descriptors are computed during this routine
# It is slow because it's currently not multi-threaded.
_AtWA, _AtWb = PotentialLearning.ooc_learn!(lb, ds_train;ws=ws,symmetrize=false, λ=0.01)
CSV.write("4body_4-8poly_fit-to-10K_POD_beta.txt", DataFrame(Tables.table(lb.β)), header=false)

# Computing the test error
println("computing test descriptors")
@time begin
edescr_test = compute_local_descriptors(ds_test, lmp_pod)
fdescr_test = compute_force_descriptors(ds_test, lmp_pod)
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
