using PotentialLearning
using CSV, DataFrames 
using Unitful, UnitfulAtomic
using InteratomicBasisPotentials

include("../../lammps_related/lammps_jl_tests/utilities/fixed_parse_extXYZ/parse_extXYZ.jl")
include("./subtract_peratom_e.jl")

function adjust_energies(ds, vref_dict)
    for config in ds
        new_energy = subtract_peratom_e(config,vref_dict)
        config.data[Energy] = new_energy
    end
end

vref_dict = Dict(:Hf => -2.70516846,
                :O => -0.01277342)

ds_path = "../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz"
ds = fixed_load_data(ds_path, ExtXYZ(u"eV", u"Å"))

train_idxs = vec(Matrix(CSV.read("Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs   = vec(Matrix(CSV.read("Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

data_train = ds[train_idxs]
data_val   = ds[val_idxs]

#data_train = ds[1:5]

ace1 = ACE(species           = [:Hf, :O],
          body_order        = 4,
          polynomial_degree = 10,
          wL                = 1.5,
          csp               = 1.0,
          r0                = 2.15,
          rcutoff           = 5.0 )

adjust_energies(data_train,vref_dict)
#adjust_energies(data_val,vref_dict)

e_descr_train = compute_local_descriptors(data_train,ace1)
f_descr_train = compute_force_descriptors(data_train,ace1)
ds_train = DataSet(data_train .+ e_descr_train .+ f_descr_train)

#lb = LBasisPotential(ace1)
lb = LBasisPotentialExt(ace1)

#learn statement 
ws, int = [1.0, 1.0], false
learn!(lb, ds_train, ws, int)

#lb, Σ = learn!(lb, ds_train; α=1e-6)

f_train = get_all_forces(ds_train)
f_train_pred = get_all_forces(ds_train,lb)
f_train_mae, f_train_rmse, f_train_rsq = calc_metrics(f_train_pred, f_train)
# (0.16270239982299103, 0.2242038939577839, 0.983185804119983)

n_train = [length(get_system(config).particles) for config in data_train]
e_train = get_all_energies(ds_train) ./ n_train
e_train_pred = get_all_energies(ds_train, lb) ./ n_train
e_mae, e_rmse, e_rsq = calc_metrics(e_train, e_train_pred)
# (0.0041530695066074215, 0.005346068866017318, 0.9994512695651514)
