using Statistics: mean 
using Distances

include("../utils/utils.jl")
includet("../utils/linear_charge_model.jl")

#param_file = "../files/4body_CH_param.pod"
#param_file = "../files/6body_CH_param.pod"
#data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
#elem_list = [:C, :H]
#target_Qtot = 0.0
#target_Qtot = 1.0
#e_lambda = 0.01

param_file = "../files/4body_AuAlMgO_param.pod"
data_fname = "../behler_charge_datasets/AuMgO/input.data"
elem_list =[:Au, :Al, :Mg, :O]
target_Qtot = 0.0
e_lambda = 0.01
#
#param_file = "../files/4body_Ag_param.pod"
#data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
#elem_list = [:Ag,]
#target_Qtot = -1.0
#target Qtot = 1.0
#e_lambda = 0.01 

#param_file = "../files/4body_NaCl_param.pod"
#data_fname = "../behler_charge_datasets/NaCl/input.data"
#elem_list = [:Na, :Cl]
#target_Qtot = 1.0
#e_lambda=1e-7


lmp_pod    = LAMMPS_POD(param_file, elem_list)
basis_size = length(lmp_pod)
@show basis_size

configs = read_behler_input(data_fname)
configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),target_Qtot, atol=1e-8)][1:1000]

test_sys = configs_train[1]
cry_to_cart = Matrix(reduce(hcat, ustrip.(bounding_box(test_sys)) )')
cart_to_cry = inv(cry_to_cart)

cart_pos = reduce(hcat, ustrip.(position(test_sys)))'
cry_pos = mapslices(x-> cart_to_cry*x, cart_pos, dims=2)

