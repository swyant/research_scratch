using CSV, DataFrames
using ACEpotentials
using JLD2 

aspirin_dsfile  = "../../../../../datasets/revMD17/rmd17_aspirin.xyz"
train_indexfile = "../../../../../datasets/revMD17/rmd17/splits/index_train_01.csv" 
test_indexfile  = "../../../../../datasets/revMD17/rmd17/splits/index_test_01.csv"

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
#train_data = data[train_indices]
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
                           totaldegree = Dict(1 => 20, 2 => 20, 3 => 20, 4 => 8),
                           wL          = 2.0,
                           transform   = (:agnesi,2,4), # default
                           envelope    = (:x, 2, 2),
                           pure2b      = false, # not default, can play around with this
                           pair_rcut   = 5.5)
                           # defaults for everything else for the pair potential

solver =  ACEfit.LSQR(damp = 0.25, atol = 1e-6)

P = smoothness_prior(combined_basis; p=1.5)

#A, Y, W = ACEfit.assemble(train_data, combined_basis)
#ddict = load("./design_matrix_v1_spl.jld2")   
#ddict = load("./mod_design_matrix_v1.jld2")
#ddict = load("./mod2_design_matrix_v1.jld2")
#A = ddict["A"]
#Ap = ddict["Ap"]
#Apw = ddict["Apw"]
#Y = ddict["Y"]
#W = ddict["W"]

#results = ACEfit.solve(solver, Apw, W .* Y)
#save("bck_results_v1.jld2", Dict("C" => results["C"]))
#save("results_v1.jld2", Dict("results" => results))

ddict = load("results_v1.jld2")
results = ddict["results"]

coeffs = P \ results["C"]

pot1 = JuLIP.MLIPs.SumIP(limited_Vref, JuLIP.MLIPs.combine(combined_basis, coeffs)) 

include("./alt_linear_errors.jl")
alt_linear_errors(test_data, pot1; energy_per_atom=false)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│ default │  10.151 │    0.036 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │  10.151 │    0.036 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│ default │   7.382 │    0.024 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   7.382 │    0.024 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
