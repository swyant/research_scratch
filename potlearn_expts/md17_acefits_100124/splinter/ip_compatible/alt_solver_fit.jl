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


#model_idx = 3
#combined_basis = ACE1x.ace_basis(elements    = elements,
#                           order       = 4,
#                           rcut        = 4.4,
#                           rin         = 0.77,
#                           r0          = my_r0,
#                           totaldegree = 9,
#                           wL          = 2.0,
#                           transform   = (:agnesi,2,4), # default
#                           envelope    = (:x, 2, 2),
#                           pure2b      = false,
#                           pair_rcut   = 5.5)

#model_idx = 5
#combined_basis = no_pair_basis(elements    = elements,
#                           order       = 4,
#                           rcut        = 4.4,
#                           rin         = 0.77,
#                           r0          = my_r0,
#                           totaldegree = 9,
#                           wL          = 2.0,
#                           transform   = (:agnesi,2,4), # default
#                           envelope    = (:x, 0, 2),
#                           pure2b      = false)
#

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

solvers = [("lsqr_damp0p25_smoothprior",
            ACEfit.LSQR(damp = 0.25, atol = 1e-6),
            smoothness_prior(combined_basis; p=1.5)),
           ("qr_lambda0p25_smoothprior",
            ACEfit.QR(lambda= 0.25),
            smoothness_prior(combined_basis; p=1.5)),
           ("qr_lambda0p01_smoothprior",
            ACEfit.QR(lambda= 0.01),
            smoothness_prior(combined_basis; p=1.5)),
           ("qr_lambda0_smoothprior",
            ACEfit.QR(lambda= 0.0),
            smoothness_prior(combined_basis; p=1.5)),
           ("qr_lambda0p25_noprior",
            ACEfit.QR(lambda= 0.25),
            nothing),
           ("qr_lambda0p01_noprior",
            ACEfit.QR(lambda= 0.01),
            nothing),
           ("qr_lambda0_noprior",
            ACEfit.QR(lambda= 0.0),
            nothing)]

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
