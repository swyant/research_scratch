using Optim
using Statistics: mean 
import Zygote
using LineSearches
includet("../utils/utils.jl")
includet("../utils/linear_charge_model.jl")

#param_file = "../files/4body_CH_param.pod"
#param_file = "../files/6body_CH_param.pod"
#data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
#elem_list = [:C, :H]
#target_Qtot = 0.0
#target_Qtot = 1.0
#e_lambda = 0.01

#param_file = "../files/4body_AuAlMgO_param.pod"
#data_fname = "../behler_charge_datasets/AuMgO/input.data"
#elem_list =[:Au, :Al, :Mg, :O]
#target_Qtot = 0.0
#e_lambda = 0.01

#param_file = "../files/4body_Ag_param.pod"
#data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
#elem_list = [:Ag,]
#target_Qtot = -1.0
#target Qtot = 1.0
#e_lambda = 0.01 

param_file = "../files/4body_NaCl_param.pod"
data_fname = "../behler_charge_datasets/NaCl/input.data"
elem_list = [:Na, :Cl]
target_Qtot = 1.0
e_lambda=1e-7


lmp_pod    = LAMMPS_POD(param_file, elem_list)
basis_size = length(lmp_pod)
@show basis_size

configs = read_behler_input(data_fname)
configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),target_Qtot, atol=1e-8)][1:1000]

# vector of vector format
qi_ref = [ustrip.(get_atomic_charges(config)) for config in configs_train]
e_ref  = [ustrip(get_energy(config)) for config in configs_train]

#ecentered_ref = e_ref .- mean(e_ref)
epa_ref = [e_ref[i]/length(configs_train[i]) for i in 1:length(configs_train)]
ecentered_ref = e_ref .- (length.(configs_train) .* mean(epa_ref))

ξ0 = learn_charge_model(configs_train, lmp_pod; λ=0.01)    

lbp_charge = LBasisChargeModel(lmp_pod,ξ0)
all_ref_charges = get_all_atomic_charges(configs_train)
all_pred_charges = get_all_atomic_charges(configs_train,lbp_charge)

charge_mae,charge_rmse,charge_mape = calc_mae_rmse_mape(all_pred_charges,all_ref_charges)
println("CHARGE\nmae: $(1000*charge_mae) me\nrmse: $(1000*charge_rmse) me\nmape: $(100*charge_mape) %\n")

α0 = learn_only_energy(configs_train, lmp_pod, ecentered_ref; λ=e_lambda)
lbp_energy = LBasisPotential(α0, [0.0,], lmp_pod)
all_pred_energies = get_all_energies(configs_train, lbp_energy) .+ (length.(configs_train) .* mean(epa_ref))
#energy_mae,energy_rmse,energy_mape = calc_mae_rmse_mape(all_pred_energies,ecentered_ref)
#println("ENERGY\nmae: $(1000*energy_mae) meV\nrmse: $(1000*energy_rmse) meV\nmape: $(100*energy_mape) %\n")

energy_mae,energy_rmse,energy_mape = calc_mae_rmse_mape(all_pred_energies ./ length.(configs_train), e_ref ./ length.(configs_train))
println("ENERGY\nmae: $(1000*energy_mae) meV per atom\nrmse: $(1000*energy_rmse) meV per atom\nmape: $(100*energy_mape) %\n")

check_lcb = LinearChargeBasis(lmp_pod,lbp_charge,3,target_Qtot,elem_list)

#params = learn_charge_energy_model(configs_train, check_lcb, ecentered_ref; λ=e_lambda)
#new_lbp = LBasisPotential(params, [0.0,], check_lcb);
#all_pred_energies2 = get_all_energies(configs_train, new_lbp) .+ (length.(configs_train) .* mean(epa_ref))
#
#energy_mae2,energy_rmse2,energy_mape2 = calc_mae_rmse_mape(all_pred_energies2 ./ length.(configs_train), e_ref ./ length.(configs_train))
#println("ENERGY\nmae: $(1000*energy_mae2) meV per atom\nrmse: $(1000*energy_rmse2) meV per atom\nmape: $(100*energy_mape2) %\n")
