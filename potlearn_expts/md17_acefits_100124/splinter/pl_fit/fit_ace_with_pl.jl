using InteratomicPotentials, PotentialLearning
using CSV, DataFrames 
using Unitful
using JLD2

include("./subtract_peratom_e.jl")

function adjust_energies(ds, vref_dict)
    for config in ds
        new_energy = subtract_peratom_e(config,vref_dict)
        config.data[Energy] = new_energy
    end
end

vref_dict = Dict(:H => -13.587222780835477,
                 :C => -1029.4889999855063,
                 :O => -2041.9816003861047)

ds_path = "../../../../../../datasets/revMD17/rmd17_aspirin.xyz"
ds = load_data(ds_path, ExtXYZ(u"eV", u"Å"))

train_indexfile = "../../../../../../datasets/revMD17/rmd17/splits/index_train_01.csv" 
test_indexfile  = "../../../../../../datasets/revMD17/rmd17/splits/index_test_01.csv"

train_idxs = vec(Matrix(CSV.read(train_indexfile,DataFrame,header=false)))
test_idxs  = vec(Matrix(CSV.read(test_indexfile,DataFrame,header=false)))

data_train  = ds[train_idxs]
data_test   = ds[test_idxs]


ace1 = ACE(species           = [:C,:O,:H],
          body_order        = 5,
          polynomial_degree = 9,
          wL                = 2.0,
          csp               = 1.0,
          r0                = 1.43,
          rcutoff           = 4.4 )

adjust_energies(data_train,vref_dict)
#adjust_energies(data_val,vref_dict)

println("computing training descriptors")
e_descr_train = compute_local_descriptors(data_train,ace1)
f_descr_train = compute_force_descriptors(data_train,ace1)
ds_train = DataSet(data_train .+ e_descr_train .+ f_descr_train)

println("computing test descriptors")
e_descr_test = compute_local_descriptors(data_test,ace1)
f_descr_test = compute_force_descriptors(data_test,ace1)
ds_test = DataSet(data_test .+ e_descr_test .+ f_descr_test)


lb = LBasisPotential(ace1)

#learn statement 
ws, int = [30/sqrt(21), 1.0], false
learn!(lb, ds_train, ws, int)

save("beta_plfit.jld2", Dict("beta" => lb.β)) 

f_test = get_all_forces(ds_test)
f_test_pred = get_all_forces(ds_test,lb)
f_test_mae, f_test_rmse, f_test_rsq = calc_metrics(f_test_pred, f_test)
@show f_test_mae, f_test_rmse
# (f_test_mae, f_test_rmse) = (0.08469358249833331, 0.24925088441993595)

e_test = get_all_energies(ds_test) 
e_test_pred = get_all_energies_w_onebody(ds_test, lb, vref_dict) 
e_test_mae, e_test_rmse, e_test_rsq = calc_metrics(e_test, e_test_pred)
@show e_test_mae, e_test_rmse
# (e_test_mae, e_test_rmse) = (0.02944178557525447, 0.05568985533122552)
