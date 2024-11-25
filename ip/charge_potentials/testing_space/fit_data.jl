
include("../utils/utils.jl")
include("../utils/charge_model.jl")

param_file = "./4body_CH_param.pod"
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

#configs = [config for config in configs 
#            if isapprox(0.0, ustrip(config[:total_charge]); atol=1e-8)]


configs = [config for config in configs 
            if isapprox(1.0, ustrip(config[:total_charge]); atol=1e-8)]

β = learn_charge_model(configs, lmp_pod)
lbp = LBasisPotential(β,[0.0,],lmp_pod)

all_ref_charges = get_all_atomic_charges(configs)
all_pred_charges = get_all_atomic_charges(configs,lbp)

mae,rmse,mape = calc_mae_rmse_mape(all_pred_charges,all_ref_charges)
println("mae: $(1000*mae) me\nrmse: $(1000*rmse) me\nmape: $(100*mape) %")
