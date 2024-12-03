using Optim

include("../utils/utils.jl")
includet("../utils/charge_model.jl")
#include("../utils/alt_charge_model.jl")

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

#param_file = "./4body_NaCl_param.pod"
#data_fname = "../behler_charge_datasets/NaCl/input.data"
#lmp_pod = LAMMPS_POD(param_file, [:Na, :Cl])

configs = read_behler_input(data_fname)
configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),-1.0)][1:1000]
qi_ref = [ustrip.(get_atomic_charges(config)) for config in configs_train]

config_set = configs_train 
βref = learn_charge_model(config_set, lmp_pod; λ=0.01)    
#βref = learn_charge_model(config_set, lmp_pod; λ=1e-10, reg_style=:scale)
lbp = LBasisPotential(βref,[0.0,],lmp_pod)

all_ref_charges = get_all_atomic_charges(config_set)
all_pred_charges = get_all_atomic_charges(config_set,lbp)

mae,rmse,mape = calc_mae_rmse_mape(all_pred_charges,all_ref_charges)
println("mae: $(1000*mae) me\nrmse: $(1000*rmse) me\nmape: $(100*mape) %\n")

ds = prepare_train_db(config_set,lmp_pod)
β0 = βref + 0.1*randn(length(lmp_pod))

obj_fun = p -> only_charge_loss(p,ds,qi_ref)

res = optimize(obj_fun, β0, LBFGS(), Optim.Options(iterations=100); autodiff=:forward)
#=
 * Status: failure (reached maximum number of iterations)

 * Candidate solution
    Final objective value:     8.340141e-02

 * Found with
    Algorithm:     L-BFGS

 * Convergence measures
    |x - x'|               = 1.10e-06 ≰ 0.0e+00
    |x - x'|/|x'|          = 3.27e-06 ≰ 0.0e+00
    |f(x) - f(x')|         = 6.75e-12 ≰ 0.0e+00
    |f(x) - f(x')|/|f(x')| = 8.10e-11 ≰ 0.0e+00
    |g(x)|                 = 5.25e-06 ≰ 1.0e-08

 * Work counters
    Seconds run:   4  (vs limit Inf)
    Iterations:    100
    f(x) calls:    362
    ∇f(x) calls:   362


=#
new_lbp = LBasisPotential(res.minimizer, [0.0,],lmp_pod)
new_pred_charges = get_all_atomic_charges(config_set,new_lbp)
@show mae,rmse,mape = calc_mae_rmse_mape(new_pred_charges,all_ref_charges)
