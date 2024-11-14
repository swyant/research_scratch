using InteratomicPotentials, PotentialLearning
using Unitful
using CSV, DataFrames
using JLD2


include("../../md17_acefits_100124/splinter/pl_fit/subtract_peratom_e.jl")

train_xyzfile = "../../../../../datasets/ISO17/my_iso17_train.xyz"

raw_train = load_data(train_xyzfile, ExtXYZ(u"eV", u"Å"))
train_idxs = vec(Matrix(CSV.read("10K_random_train_idxs",DataFrame,header=false)))

ds_train = raw_train[train_idxs]
save("10K_train_configs.jld2", Dict("ds_train" => ds_train))
