using Optim
using Statistics: mean 
import Zygote
using LineSearches

includet("../utils/utils.jl")
includet("../utils/linear_charge_model.jl")

#param_file = "../files/4body_CH_param.pod"
##param_file = "../files/6body_CH_param.pod"
#data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
#elem_list = [:C, :H]
##target_Qtot = 0.0
#target_Qtot = 1.0
#e_lambda = 0.01
#rcut = 10.0

param_file = "../files/4body_AuAlMgO_param.pod"
data_fname = "../behler_charge_datasets/AuMgO/input.data"
elem_list =[:Au, :Al, :Mg, :O]
target_Qtot = 0.0
e_lambda = 0.01
rcut = 4.5

#param_file = "../files/4body_Ag_param.pod"
#data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
#elem_list = [:Ag,]
##target_Qtot = -1.0
#target_Qtot = 1.0
#e_lambda = 0.01 
#rcut = 10.0

#param_file = "../files/4body_NaCl_param.pod"
#data_fname = "../behler_charge_datasets/NaCl/input.data"
#elem_list = [:Na, :Cl]
#target_Qtot = 1.0
#e_lambda=1e-7
#rcut = 10.0

lmp_pod    = LAMMPS_POD(param_file, elem_list)
basis_size = length(lmp_pod)
@show basis_size

configs = read_behler_input(data_fname)
configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),target_Qtot, atol=1e-8)][1:1000]

ξ0 = learn_charge_model(configs_train, lmp_pod; λ=0.01)    

lbp_charge = LBasisChargeModel(lmp_pod,ξ0)
all_ref_charges = get_all_atomic_charges(configs_train)
all_pred_charges = get_all_atomic_charges(configs_train,lbp_charge)

charge_mae,charge_rmse,charge_mape = calc_mae_rmse_mape(all_pred_charges,all_ref_charges)
println("CHARGE\nmae: $(1000*charge_mae) me\nrmse: $(1000*charge_rmse) me\nmape: $(100*charge_mape) %\n")

 
test_sys = configs_train[1]
ref_charge_derivs = finite_difference_charges(test_sys,lbp_charge)
charge_derivs = compute_atomic_charge_derivatives(test_sys,lbp_charge)

