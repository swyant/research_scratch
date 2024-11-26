
include("../utils/utils.jl")
include("../utils/charge_model.jl")
#include("../utils/alt_charge_model.jl")

param_file = "./4body_CH_param.pod"
#param_file = "./6body_CH_param.pod"
data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
lmp_pod = LAMMPS_POD(param_file, [:C, :H])

#param_file = "./4body_AuAlMgO_param.pod"
#data_fname = "../behler_charge_datasets/AuMgO/input.data"
#lmp_pod = LAMMPS_POD(param_file, [:Au, :Al, :Mg, :O])

#param_file = "./4body_Ag_param.pod"
#data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
#lmp_pod    = LAMMPS_POD(param_file, [:Ag,])

#param_file = "./4body_NaCl_param.pod"
#data_fname = "../behler_charge_datasets/NaCl/input.data"
#lmp_pod = LAMMPS_POD(param_file, [:Na, :Cl])

configs = read_behler_input(data_fname)

configs0 = [config for config in configs 
            if isapprox(0.0, ustrip(config[:total_charge]); atol=1e-8)]


configs1 = [config for config in configs 
            if isapprox(1.0, ustrip(config[:total_charge]); atol=1e-8)]

config_dict = Dict("all_data" => configs,
                   "0 charge" => configs0,
                   "+1 charge"=> configs1)

for key in keys(config_dict)
    config_set = config_dict[key]
    #β = learn_charge_model(config_set, lmp_pod; λ=0.001)
    
    β = learn_charge_model(config_set, lmp_pod; λ=1e-10, reg_style=:scale)
    lbp = LBasisPotential(β,[0.0,],lmp_pod)
    
    all_ref_charges = get_all_atomic_charges(config_set)
    all_pred_charges = get_all_atomic_charges(config_set,lbp)
    
    mae,rmse,mape = calc_mae_rmse_mape(all_pred_charges,all_ref_charges)
    println("configs set: $(key)")
    println("mae: $(1000*mae) me\nrmse: $(1000*rmse) me\nmape: $(100*mape) %\n")
end
