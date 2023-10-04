
using Distributed; addprocs(4)

using ACE1, ACE1pack, ACEfit
using JuLIP
using CSV, DataFrames 
using JLD2
using Random 

include("../../utilities/custom_export2lammps/custom_export2lammps.jl")

function obtain_train_idxs(train_data,frac)
    @assert frac <= 1.0
    total_num_data = length(train_data)
    num_train = Int(floor(frac*total_num_data))

    perm_idxs = Random.randperm(total_num_data)
    train_idxs = perm_idxs[begin:1:num_train]
    train_idxs
end

ds_path = "../../../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz"
raw_data = read_extxyz(ds_path)
label = "HfO2_N2_P6_r4"

rpib = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 2, 
            maxdeg  = 6,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 4.0,
            pin     = 0,
)

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))
vref = JuLIP.OneBody([:Hf => -2.70516846, :O => -0.01277342]...)

all_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data]
#=
- for each ace 
- select random indices, save those indices to a parameterized file 
- make dataset for those random indices (can I do it from all_data?)
=#

all_train_idxs = vec(Matrix(CSV.read("./indices/Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs = vec(Matrix(CSV.read("./indices/Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

fit_idx = 1
train_idxs = obtain_train_idxs(all_data[all_train_idxs], 0.9)
CSV.write("indices/fit_$(fit_idx)_train_idxs.csv",DataFrame(Tables.table(train_idxs)),header=false)

A, Y, W = ACEfit.assemble(all_data[train_idxs], rpib)
solver_reg = ACEfit.QR(; lambda=1e-3)
results_reg = ACEfit.solve(solver_reg,A,Y)
CSV.write("fits/$(label)_fit$(fit_idx)_coeffs.csv", DataFrame(Tables.table(results_reg["C"])),header=false)

pot_reg_tmp = JuLIP.MLIPs.combine(rpib,results_reg["C"])
pot_reg = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
train_errors = ACE1pack.linear_errors(all_data[train_idxs],pot_reg)
val_errors = ACE1pack.linear_errors(all_data[val_idxs],pot_reg)

save("./errors/$(label)_fit$(fit_idx)_errors.jld2","train_errors",train_errors,"val_errors",val_errors)
custom_export2lammps("./fits/$(label)_fit$(fit_idx).yace", pot_reg,rpib)
