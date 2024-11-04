using InteratomicPotentials, PotentialLearning
using CSV, DataFrames
using Unitful
using JLD2 
using LinearAlgebra 

include("./no_pair_basis.jl") # also using ACEpotentials

aspirin_xyzfile = "../../../../../../datasets/revMD17/rmd17_aspirin.xyz"

train_indexfile = "../../../../../../datasets/revMD17/rmd17/splits/index_train_01.csv" 
test_indexfile  = "../../../../../../datasets/revMD17/rmd17/splits/index_test_01.csv"

train_idxs = vec(Matrix(CSV.read(train_indexfile,DataFrame,header=false)))
test_idxs  = vec(Matrix(CSV.read(test_indexfile,DataFrame,header=false)))


ace = ACE(species           = [:C,:O,:H],
      body_order        = 3,
      polynomial_degree = 10,
      wL                = 2.0,
      csp               = 1.0,
      r0                = 1.43,
      rcutoff           = 4.4 )

@show length(ace)

include("../../splinter/pl_fit/subtract_peratom_e.jl")

lb = LBasisPotential(ace)

vref_dict = Dict(:H => -13.587222780835477,
                 :C => -1029.4889999855063,
                 :O => -2041.9816003861047)


ds = load_data(aspirin_xyzfile, ExtXYZ(u"eV", u"Å"))
data_train  = ds[train_idxs]
data_test   = ds[test_idxs]

adjust_energies(data_train,vref_dict)

println("computing training descriptors")
e_descr_train = compute_local_descriptors(data_train,ace)
f_descr_train = compute_force_descriptors(data_train,ace)
ds_train = DataSet(data_train .+ e_descr_train .+ f_descr_train)

println("computing test descriptors")
e_descr_test = compute_local_descriptors(data_test,ace)
f_descr_test = compute_force_descriptors(data_test,ace)
ds_test = DataSet(data_test .+ e_descr_test .+ f_descr_test)

e_test = get_all_energies(ds_test)
f_test = get_all_forces(ds_test)


#ws, int = [30.0, 1.0], false
ws, int = [1.0, 1.0], false
lp = PotentialLearning.LinearProblem(ds_train)


#A, W, b = PotentialLearning.assemble_matrices(lp,ws)
#learn!(lp,ws,int)
#default_learn_β = lp.β[:]
#
#original_AtWA, original_AtWb  = PotentialLearning.ooc_learn!(lb, ds_train; ws=ws, symmetrize=false, lambda=nothing)
#
#vanilla_ooc_β = lb.β[:]
#
#_AtWA, _AtWb = PotentialLearning.ooc_learn!(lb, ds_train; 
#                                            ws=ws, 
#                                            reg_style=:default, 
#                                            symmetrize=false, 
#                                            lambda=0.01,
#                                            AtWA = original_AtWA,
#                                            AtWb = original_AtWb)
#
#default_reg_β = lb.β[:]
#
#_AtWA, _AtWb = PotentialLearning.ooc_learn!(lb, ds_train; 
#                                            ws=ws, 
#                                            reg_style=:scale_thresh, 
#                                            symmetrize=true, 
#                                            lambda=1e-10,
#                                            AtWA = original_AtWA,
#                                            AtWb = original_AtWb)
#
#cuong_reg_β = lb.β[:]
#
#save("beta_sets.jld2", Dict("default_learn_β" => default_learn_β,
#                            "vanilla_ooc_β"   => vanilla_ooc_β,
#                            "default_reg_β"   => default_reg_β,
#                            "cuong_reg_β"     => cuong_reg_β))
