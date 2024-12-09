using Optim
using Statistics: mean 

include("../utils/utils.jl")
includet("../utils/charge_model.jl")

#param_file = "./4body_CH_param.pod"
#param_file = "./6body_CH_param.pod"
#data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
#lmp_pod = LAMMPS_POD(param_file, [:C, :H])

#param_file = "./4body_AuAlMgO_param.pod"
#data_fname = "../behler_charge_datasets/AuMgO/input.data"
#lmp_pod = LAMMPS_POD(param_file, [:Au, :Al, :Mg, :O])

param_file = "./4body_Ag_param.pod"
data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
lmp_pod    = LAMMPS_POD(param_file, [:Ag,])
basis_size = length(lmp_pod)

#param_file = "./4body_NaCl_param.pod"
#data_fname = "../behler_charge_datasets/NaCl/input.data"
#lmp_pod = LAMMPS_POD(param_file, [:Na, :Cl])

configs = read_behler_input(data_fname)
configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),-1.0)][1:1000]

# vector of vector format
qi_ref = [ustrip.(get_atomic_charges(config)) for config in configs_train]
e_ref  = [ustrip(get_energy(config)) for config in configs_train]
ecentered_ref = e_ref .- mean(e_ref)

β0 = learn_charge_model(configs_train, lmp_pod; λ=0.01)    

lbp_charge = LBasisPotential(β0,[0.0,],lmp_pod)
all_ref_charges = get_all_atomic_charges(configs_train)
all_pred_charges = get_all_atomic_charges(configs_train,lbp_charge)

charge_mae,charge_rmse,charge_mape = calc_mae_rmse_mape(all_pred_charges,all_ref_charges)
println("CHARGE\nmae: $(1000*charge_mae) me\nrmse: $(1000*charge_rmse) me\nmape: $(100*charge_mape) %\n")

α0 = learn_only_energy(configs_train, lmp_pod, ecentered_ref; λ=0.01)
lbp_energy = LBasisPotential(α0, [0.0,], lmp_pod)
all_pred_energies = get_all_energies(configs_train, lbp_energy)

energy_mae,energy_rmse,energy_mape = calc_mae_rmse_mape(all_pred_energies,ecentered_ref)
println("ENERGY\nmae: $(1000*energy_mae) meV\nrmse: $(1000*energy_rmse) meV\nmape: $(100*energy_mape) %\n")

γ0 = randn(basis_size)

#ps0 = [α0+randn(basis_size);β0+randn(basis_size);γ0]
ps0 = [β0+randn(basis_size);γ0]

cbp0 = ChargedBasisPotential(lmp_pod, zeros(basis_size), zeros(basis_size), zeros(basis_size))

ds = prepare_train_db(configs_train,lmp_pod)
obj_fun = p -> simpler_loss(p,ds, cbp0, qi_ref, ecentered_ref; we=30.0, λ=0.01)

#res = optimize(obj_fun, ps0, LBFGS(), Optim.Options(iterations=1000); autodiff=:forward)
