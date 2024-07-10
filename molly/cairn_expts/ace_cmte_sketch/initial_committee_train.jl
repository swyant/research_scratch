using PotentialLearning, InteratomicPotentials 
using Unitful 
using CSV, DataFrames

cd("./pce_sketch_2")
include("./training_utils.jl")

# ace specification 
ace = ACE(  species            = [:Hf, :O],
            body_order         = 4, 
            polynomial_degree  = 10,
            wL                 = 1.5, 
            csp                = 1.0,
            r0                 = 2.15,
            rcutoff            = 5.0)

# load up data sets
ds_path  = "/Users/swyant/cesmix/datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz"
all_train_idxs = vec(Matrix(CSV.read("./Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))

all_confs = load_data(ds_path, ExtXYZ(u"eV", u"Å" ))
train_confs = all_confs[all_train_idxs]
e_descr_train = compute_local_descriptors(train_confs,ace)
f_descr_train = compute_force_descriptors(train_confs,ace)

ds_train = DataSet(train_confs .+ e_descr_train .+ f_descr_train)

lbps = perform_and_record_simple_fits(ds_train,ace,10)

#train_idxs = obtain_train_idxs(0.6,length(ds_train))
#lbp_trial = LBasisPotential(ace)
#ws, int = [1000.0, 1.0], false
#learn!(lbp_trial,ds_train[train_idxs],ws,int)
#
#CSV.write("trial_coeffs.csv", DataFrame(Tables.table(lbp_trial.β)),header=false)
 