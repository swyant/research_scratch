using InteratomicPotentials, PotentialLearning
using CSV, DataFrames
using Unitful
using JLD2 
using LinearAlgebra: norm
using Statistics: mean

aspirin_xyzfile = "../../../../../../datasets/revMD17/rmd17_aspirin.xyz"

train_indexfile = "../../../../../../datasets/revMD17/rmd17/splits/index_train_01.csv" 
test_indexfile  = "../../../../../../datasets/revMD17/rmd17/splits/index_test_01.csv"

train_idxs = vec(Matrix(CSV.read(train_indexfile,DataFrame,header=false)))
test_idxs  = vec(Matrix(CSV.read(test_indexfile,DataFrame,header=false)))

ace = ACE(species           = [:C,:O,:H],
      body_order        = 5,
      polynomial_degree = 9,
      wL                = 2.0,
      csp               = 1.0,
      r0                = 1.43,
      rcutoff           = 4.4 )

#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 5,
#      polynomial_degree = 8,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )

#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 5,
#      polynomial_degree = 7,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )


#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 5,
#      polynomial_degree = 6,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )

#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 4,
#      polynomial_degree = 9,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )

#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 4,
#      polynomial_degree = 8,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )

#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 3,
#      polynomial_degree = 9,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )

#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 3,
#      polynomial_degree = 12,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )


#ace = ACE(species           = [:C,:O,:H],
#      body_order        = 3,
#      polynomial_degree = 10,
#      wL                = 2.0,
#      csp               = 1.0,
#      r0                = 1.43,
#      rcutoff           = 4.4 )

include("../../splinter/pl_fit/subtract_peratom_e.jl")

lb = LBasisPotential(ace)

ds = load_data(aspirin_xyzfile, ExtXYZ(u"eV", u"Å"))
data_train  = ds[train_idxs]
data_test   = ds[test_idxs]

e_train = get_all_energies(data_train)
num_atoms = length(get_system(data_train[1])) # all systems have same number of atoms, md17-aspirin

avg_energy_per_atom = mean(e_train)/num_atoms
vref_dict = Dict(:H => avg_energy_per_atom,
                 :C => avg_energy_per_atom,
                 :O => avg_energy_per_atom)

adjust_energies(data_train,vref_dict)

println("computing training descriptors")
e_descr_train = compute_local_descriptors(data_train,ace)
f_descr_train = compute_force_descriptors(data_train,ace)
ds_train = DataSet(data_train .+ e_descr_train .+ f_descr_train)

println("computing test descriptors")
e_descr_test = compute_local_descriptors(data_test,ace)
f_descr_test = compute_force_descriptors(data_test,ace)
ds_test = DataSet(data_test .+ e_descr_test .+ f_descr_test)

@show length(ace)

ws, int = [30.0, 1.0], true
lp = PotentialLearning.LinearProblem(ds_train)
learn!(lp,ws,int; λ=0.01)
lb.β .= lp.β
lb.β0 .= lp.β0

e_test = get_all_energies(ds_test)
f_test = get_all_forces(ds_test)

###### evaluate test errors 
#f_test_pred = get_all_forces(ds_test,lb)
#f_test_mae, f_test_rmse, f_test_rsq = calc_metrics(f_test_pred, f_test)
#@show f_test_mae, f_test_rmse
#
#e_test_pred = get_all_energies_w_onebody(ds_test, lb, vref_dict)
#e_test_mae, e_test_rmse, e_test_rsq = calc_metrics(e_test, e_test_pred)
#@show e_test_mae, e_test_rmse
#@show norm(lb.β)
