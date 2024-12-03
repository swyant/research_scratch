
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

configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),-1.0)][1:1000]
#configs_train = [config for config in configs if isapprox(ustrip(get_total_charge(config)),-1.0)][1:1000]

qi_ref = [ustrip.(get_atomic_charges(config)) for config in config_train]
e_ref = [ustrip(get_energy(config)) for config in config_train]

e_train = get_all_energies(base_ds_train)
e_test  = get_all_energies(base_ds_test)
f_test  = get_all_forces(base_ds_test)

num_atoms = length(get_system(base_ds_train[1])) # all systems have same number of atoms, entire ISO17 database

avg_energy_per_atom = mean(vcat(e_train,e_test))/num_atoms
vref_dict = Dict(:H => avg_energy_per_atom,
                 :C => avg_energy_per_atom,
                 :O => avg_energy_per_atom)

adjust_energies(base_ds_train,vref_dict)



