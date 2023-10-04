
using Distributed; addprocs(4)

using ACE1, ACE1pack, ACEfit
using JuLIP
using CSV, DataFrames 
using JLD2
using Random 

include("../../utilities/custom_export2lammps/custom_export2lammps.jl")

#= 
Why do it this way: I may want to naively take a random subset of the full data, or 
 I already have a list of subset idx, and I'm taking an additional subset of that, but 
 I still want the idx to be wrt the original dataset
 =#
function obtain_train_idxs(frac; train_subset=none, set_size=nothing)
    @assert frac <= 1.0
    if !isnothing(train_subset)
        if !isnothing(set_size)
            print("set size will be overwritten")
        end
        set_size = length(train_subset)
    end
    num_select = Int(floor(frac*set_size))

    perm_idxs = Random.randperm(set_size)
    rand_set_idxs = perm_idxs[begin:1:num_select]

    if isnothing(train_subset)
        return rand_set_idxs
    else 
        return train_subset[rand_set_idxs]
    end
end

function perform_and_record_fits(full_ds, rpib, vref, label, num_fits; frac=0.9, train_subset=nothing, val_subset=nothing)
    for fit_idx in 1:num_fits
        print("###### RUNNING FIT $(fit_idx) ######")
        if isnothing(train_subset)
            train_idxs = obtain_train_idxs(frac; set_size=length(full_ds))
        else
            train_idxs = obtain_train_idxs(frac; train_subset=train_subset)
        end

        CSV.write("indices/fit_$(fit_idx)_train_idxs.csv",DataFrame(Tables.table(train_idxs)),header=false)
    
        A, Y, W = ACEfit.assemble(full_ds[train_idxs], rpib)
        solver_reg = ACEfit.QR(; lambda=1e-3)
        results_reg = ACEfit.solve(solver_reg,A,Y)
        CSV.write("fits/$(label)_fit$(fit_idx)_coeffs.csv", DataFrame(Tables.table(results_reg["C"])),header=false)
    
        pot_reg_tmp = JuLIP.MLIPs.combine(rpib,results_reg["C"])
        pot_reg = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
        train_errors = ACE1pack.linear_errors(full_ds[train_idxs],pot_reg)
        if !isnothing(val_subset)
            val_errors = ACE1pack.linear_errors(full_ds[val_subset],pot_reg)
        end
    
        save("./errors/$(label)_fit$(fit_idx)_errors.jld2","train_errors",train_errors,"val_errors",val_errors)
        custom_export2lammps("./fits/$(label)_fit$(fit_idx).yace", pot_reg,rpib)
    end
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

all_train_idxs = vec(Matrix(CSV.read("./indices/Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs = vec(Matrix(CSV.read("./indices/Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

perform_and_record_fits(all_data,rpib,vref,label,5;train_subset=all_train_idxs, val_subset=val_idxs)