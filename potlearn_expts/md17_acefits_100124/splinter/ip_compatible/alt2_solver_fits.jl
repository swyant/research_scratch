using CSV, DataFrames
using ACEpotentials
using JLD2 
using LinearAlgebra

aspirin_dsfile  = "../../../../../../datasets/revMD17/rmd17_aspirin.xyz"
train_indexfile = "../../../../../../datasets/revMD17/rmd17/splits/index_train_01.csv" 
test_indexfile  = "../../../../../../datasets/revMD17/rmd17/splits/index_test_01.csv"

train_indices = vec(Matrix(CSV.read(train_indexfile,DataFrame,header=false)))
test_indices  = vec(Matrix(CSV.read(test_indexfile,DataFrame,header=false)))

raw_data = read_extxyz(aspirin_dsfile)

weights = Dict("default" => Dict("E" => 30.0, "F" => 1.0, "V" => 1.0))

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

include("./no_pair_basis.jl")

my_r0 = 1.43


# Number 6
model_idx = 6
combined_basis = no_pair_basis(elements    = elements,
                           order       = 4,
                           rcut        = 4.4,
                           rin         = 0.77,
                           r0          = my_r0,
                           totaldegree = 9,
                           wL          = 2.0,
                           transform   = PolyTransform(2,my_r0), # default
                           envelope    = (:x, 0, 2),
                           pure2b      = false)


ddict = load("design_matrix$(model_idx).jld2")
A = ddict["A"]
Y = ddict["Y"]
W = ddict["W"]

solvers = [#("lsqr_damp0p25_smoothprior",
           # ACEfit.LSQR(damp = 0.25, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p1_smoothprior",
           # ACEfit.LSQR(damp = 0.1, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p05_smoothprior",
           # ACEfit.LSQR(damp = 0.05, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p01_smoothprior",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p005_smoothprior",
           # ACEfit.LSQR(damp = 0.005, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p001_smoothprior",
           # ACEfit.LSQR(damp = 0.001, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p0005_smoothprior",
           # ACEfit.LSQR(damp = 0.0005, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p0001_smoothprior",
           # ACEfit.LSQR(damp = 0.0001, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p01_smoothprior_0p75",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=0.75)),            
           # ("lsqr_damp0p01_smoothprior_1p0",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.0)), 
           # ("lsqr_damp0p01_smoothprior_1p25",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.25)),  
           # ("lsqr_damp0p01_smoothprior_1p5",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)), 
           # ("lsqr_damp0p01_smoothprior_1p75",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.75)), 
           # ("lsqr_damp0p01_smoothprior_2p0",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=2.0)), 
           # ("lsqr_damp0p01_smoothprior_2p25",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # smoothness_prior(combined_basis; p=2.25)), 
           # ("lsqr_damp0_smoothprior",
           # ACEfit.LSQR(damp = 0.0, atol = 1e-6),
           # smoothness_prior(combined_basis; p=1.5)),
           # ("lsqr_damp0p01_noprior",
           # ACEfit.LSQR(damp = 0.01, atol = 1e-6),
           # nothing),
           # ("lsqr_damp0_noprior",
           # ACEfit.LSQR(damp = 0.0, atol = 1e-6),
           # nothing),
            ("qr_lambda0_smoothprior",
            ACEfit.QR(lambda= 0.0),
            smoothness_prior(combined_basis; p=1.5)),          
          ]

for (solver_id, solver, P) in solvers
    if !isnothing(P)
        Apw = Diagonal(W) * (A / P)
    else
        Apw = Diagonal(W)*A
    end
    
    results = ACEfit.solve(solver, Apw, W .* Y)
    
    save("results$(model_idx)_solver_$(solver_id).jld2", Dict("results" => results))
    
    if !isnothing(P)
        coeffs = P \ results["C"]
    else
        coeffs = results["C"]
    end
    
    pot = JuLIP.MLIPs.SumIP(limited_Vref, JuLIP.MLIPs.combine(combined_basis, coeffs)) 
    
    include("../alt_linear_errors.jl")

    println("SOLVER: $(solver_id)")
    errors = alt_linear_errors(test_data, pot; energy_per_atom=false)
    
    save("errors$(model_idx)_solver_$(solver_id).jld2", Dict("errors" => errors))
end
