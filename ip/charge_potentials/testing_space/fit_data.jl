
include("../utils/utils.jl")
include("../utils/charge_model.jl")

param_file = "./4body_CH_param.pod"

#data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
#data_fname = "../behler_charge_datasets/AuMgO/input.data"
data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
#data_fname = "../behler_charge_datasets/NaCl/input.data"

configs = read_behler_input(data_fname)

zero_charge_configs = [config for config in configs 
                       if isapprox(0.0, ustrip(config[:total_charge]); atol=1e-8)]

lmp_pod = LAMMPS_POD(param_file, [:C, :H])

β = learn_charge_model(zero_charge_configs, lmp_pod)
lbp = LBasisPotential(β,[0.0,],lmp_pod)

all_ref_charges = get_all_atomic_charges(zero_charge_configs)
all_pred_charges = get_all_atomic_charges(zero_charge_configs,lbp)

mae,rmse,mape = calc_mae_rmse_mape(all_pred_charges,all_ref_charges)
println("mae: $(mae) e\nrmse: $(rmse) e\nmape: $(100*mape) %")
