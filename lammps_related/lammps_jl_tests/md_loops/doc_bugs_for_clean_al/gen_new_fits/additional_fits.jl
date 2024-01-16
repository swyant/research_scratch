using Distributed; addprocs(4)

# load things 
include("./new_code/dependencies.jl")
include("./new_code/md_active_learning.jl")
include("./new_code/train_utils.jl")

# initial parameters 
ds_path = "../../../../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz"
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

initial_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data]

all_train_idxs = vec(Matrix(CSV.read("./Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs = vec(Matrix(CSV.read("./Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

perform_and_record_fits(initial_data,rpib,vref,label,15;train_subset=all_train_idxs, val_subset=val_idxs)
