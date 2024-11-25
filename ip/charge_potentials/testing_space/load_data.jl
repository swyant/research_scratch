
include("../utils/utils.jl")
#data_fname = "../behler_charge_datasets/Ag_cluster/input.data"
#data_fname = "../behler_charge_datasets/AuMgO/input.data"
data_fname = "../behler_charge_datasets/Carbon_chain/input.data"
#data_fname = "../behler_charge_datasets/NaCl/input.data"

configs = read_behler_input(data_fname)
