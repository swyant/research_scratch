using CSV, DataFrames
using ACEpotentials
using JLD2 

aspirin_dsfile  = "../../../../../../datasets/revMD17/rmd17_aspirin.xyz"
train_indexfile = "../../../../../../datasets/revMD17/rmd17/splits/index_train_01.csv" 
test_indexfile  = "../../../../../../datasets/revMD17/rmd17/splits/index_test_01.csv"

train_indices = vec(Matrix(CSV.read(train_indexfile,DataFrame,header=false)))
test_indices  = vec(Matrix(CSV.read(test_indexfile,DataFrame,header=false)))

raw_data = read_extxyz(aspirin_dsfile)

weights = Dict("default" => Dict("E" => 30.0, "F" => 1.0, "V" => 1.0))

# From https://github.com/ACEsuit/mace-jax/blob/4b899de2101c6e2085ee972aeac0e46a334fd9a0/configs/aspirin.gin#L38
Vref = OneBody(:H => -13.587222780835477,
               :C => -1029.4889999855063,
               :N => -1484.9814568572233,
               :O => -2041.9816003861047)

data = [ AtomsData(at; energy_key = "energy", force_key = "forces", 
                   weights=weights, v_ref=Vref) for at in raw_data]
train_data = data[train_indices]
test_data  = data[test_indices]

#
elements = [:C,:O,:H]


limited_Vref = OneBody(:H => -13.587222780835477,
                       :C => -1029.4889999855063,
                       :O => -2041.9816003861047)

combined_basis = ACE1x.ace_basis(elements    = elements,
                           order       = 4,
                           rcut        = 4.4,
                           rin         = 0.77,
                           r0          = :bondlen,
                           totaldegree = Dict(1 => 14, 2 => 14, 3 => 14, 4 => 8),
                           wL          = 2.0,
                           transform   = (:agnesi,2,4), # default
                           envelope    = (:x, 2, 2),
                           pure2b      = false, # not default, can play around with this
                           pair_rcut   = 5.5)
                           # defaults for everything else for the pair potential

solver =  ACEfit.LSQR(damp = 0.25, atol = 1e-6)

P = smoothness_prior(combined_basis; p=1.5)

A, Y, W = ACEfit.assemble(train_data, combined_basis)

save("design_matrix_v2.jld2", Dict("A" => A, "Y" => Y, "W" => W))
