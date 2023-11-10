using Distributed; addprocs(4)

include("./new_code/"*
        "dependencies.jl")
include("./new_code/train_utils.jl")


ds_train_path = "./pod-hfo2/train/HfO2_figshare_form_sorted_train.extxyz"

train_data = read_extxyz(ds_train_path)

rpib = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 6, 
            maxdeg  = 12,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 5.5,
            pin     = 0,
)

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))
vref = JuLIP.OneBody([:Hf => -2.70516846, :O => -0.01277342]...)


initial_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in train_data]

A, Y, W = ACEfit.assemble(initial_data, rpib)
solver_reg = ACEfit.QR(; lambda=1e-3)
results_reg = ACEfit.solve(solver_reg,A,Y)