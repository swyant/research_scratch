
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

basis_size = length(lmp_pod)
test_cbp = ChargedBasisPotential(lmp_pod, 10*randn(basis_size), 10*randn(basis_size), 10*randn(basis_size))

test_config = configs[1]
test_qtot = ustrip(get_total_charge(test_config))
