using InteratomicPotentials, PotentialLearning
using Unitful
using CSV, DataFrames
using Random
using JLD2
using LinearAlgebra: norm
using Statistics: mean 

include("../../md17_acefits_100124/splinter/pl_fit/subtract_peratom_e.jl")

train_xyzfile = "../../../../../datasets/ISO17/my_iso17_train.xyz"
test_xyzfile  = "../../../../../datasets/ISO17/my_iso17_test.xyz"

#basis_idx = 0 
basis_idx = 1

### 0
#ace = ACE(species           = [:C,:O,:H],
#          body_order        = 5,
#          polynomial_degree = 9,
#          wL                = 2.0,
#          csp               = 1.0,
#          r0                = 1.43,
#          rcutoff           = 4.4 )

## 1
ace = ACE(species           = [:C,:O,:H],
          body_order        = 3,
          polynomial_degree = 12,
          wL                = 2.0,
          csp               = 1.0,
          r0                = 1.43,
          rcutoff           = 4.4 )

@show length(ace)

raw_train = load_data(train_xyzfile, ExtXYZ(u"eV", u"Å"))
raw_test  = load_data(test_xyzfile, ExtXYZ(u"eV", u"Å"))

#train_idxs = Random.randperm(length(raw_train))[1:10000]
#test_idxs  = Random.randperm(length(raw_test))[1:10000]
#CSV.write("10K_random_train_idxs", DataFrame(Tables.table(train_idxs)),header=false)
#CSV.write("10K_random_test_idxs", DataFrame(Tables.table(test_idxs)),header=false)

train_idxs = vec(Matrix(CSV.read("10K_random_train_idxs",DataFrame,header=false)))
#test_idxs  = vec(Matrix(CSV.read("10K_random_test_idxs",DataFrame,header=false)))[1:500]
test_idxs  = vec(Matrix(CSV.read("10K_random_test_idxs",DataFrame,header=false)))

ds_train = raw_train[train_idxs]
ds_test  = raw_test[test_idxs]

e_train = get_all_energies(ds_train)
e_test = get_all_energies(ds_test)

num_atoms = length(get_system(ds_train[1])) # all systems have same number of atoms, entire ISO17 database

#avg_energy_per_atom = mean(e_train)/num_atoms
avg_energy_per_atom = mean(vcat(e_train,e_test))/num_atoms
vref_dict = Dict(:H => avg_energy_per_atom,
                 :C => avg_energy_per_atom,
                 :O => avg_energy_per_atom)

adjust_energies(ds_train,vref_dict)


#edescr_train = compute_local_descriptors(ds_train, ace)
#fdescr_train = compute_force_descriptors(ds_train, ace)
#save("train_local_descriptors$(basis_idx).jld2", Dict("edescr_train" => edescr_train))
#save("train_force_descriptors$(basis_idx).jld2", Dict("fdescr_train" =>fdescr_train))

#print("computing test descriptors")
#edescr_test = compute_local_descriptors(ds_test, ace)
#fdescr_test = compute_force_descriptors(ds_test, ace)
#save("test_local_descriptors$(basis_idx).jld2", Dict("edescr_test" => edescr_test))
#save("test_force_descriptors$(basis_idx).jld2", Dict("fdescr_test" => fdescr_test))

edescr_train = load("train_local_descriptors$(basis_idx).jld2", "edescr_train")
fdescr_train = load("train_force_descriptors$(basis_idx).jld2", "fdescr_train")
ds_train = DataSet(ds_train .+ edescr_train .+ fdescr_train)

edescr_test = load("test_local_descriptors$(basis_idx).jld2", "edescr_test")
fdescr_test = load("test_force_descriptors$(basis_idx).jld2", "fdescr_test")
ds_test = DataSet(ds_test .+ edescr_test .+ fdescr_test)

lp = PotentialLearning.LinearProblem(ds_train)
lb = LBasisPotential(ace)
###### fitting procedure
println("fitting")
#ws, int = [30.0, 1.0], true
ws, int = [30.0, 1.0], false
learn!(lp,ws,int; λ=0.01)
lb.β .= lp.β
#lb.β0 .= lp.β0


#print("computing test descriptors")
#edescr_test = compute_local_descriptors(ds_test, ace)
#fdescr_test = compute_force_descriptors(ds_test, ace)
#ds_test = DataSet(ds_test .+ edescr_test .+ fdescr_test)


f_test = get_all_forces(ds_test)
###### evaluate test errors 
#f_test_pred = get_all_forces(ds_test,lb)
#f_test_mae, f_test_rmse, f_test_rsq = calc_metrics(f_test_pred, f_test)
#@show f_test_mae, f_test_rmse
#
#e_test_pred = get_all_energies(ds_test, lb)
#e_test_mae, e_test_rmse, e_test_rsq = calc_metrics(e_test, e_test_pred)
#@show e_test_mae, e_test_rmse
#@show norm(lb.β)
