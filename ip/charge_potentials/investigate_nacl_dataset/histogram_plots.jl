using Plots
using Statistics: mean 

include("../utils/utils.jl")
includet("../utils/linear_charge_model.jl")

param_file = "../files/4body_CH_param.pod"
param_file = "../files/6body_CH_param.pod"
data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
elem_list = [:C, :H]
target_Qtot = 0.0
#target_Qtot = 1.0

#param_file = "../files/4body_AuAlMgO_param.pod"
#data_fname = "../behler_charge_datasets/AuMgO/input.data"
#elem_list =[:Au, :Al, :Mg, :O]
#target_Qtot = 0.0

#param_file = "../files/4body_Ag_param.pod"
#data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
#elem_list = [:Ag,]
#target_Qtot = -1.0
#target Qtot = 1.0

#param_file = "../files/4body_NaCl_param.pod"
#data_fname = "../behler_charge_datasets/NaCl/input.data"
#elem_list = [:Na, :Cl]
#target_Qtot = 1.0

configs = read_behler_input(data_fname)
configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),target_Qtot, atol=1e-8)][1:1000]

# vector of vector format
qi_ref = [ustrip.(get_atomic_charges(config)) for config in configs_train]
e_ref  = [ustrip(get_energy(config)) for config in configs_train]

epa_ref = [e_ref[i]/length(configs_train[i]) for i in 1:length(configs_train)]
ecentered_ref = e_ref .- (length.(configs_train) .* mean(epa_ref))

