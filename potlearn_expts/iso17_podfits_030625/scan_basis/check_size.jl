using InteratomicPotentials, PotentialLearning
using Unitful
using Random
using JLD2
using LinearAlgebra: norm
using Statistics: mean 

include("../files/subtract_peratom_e.jl")

function compute_error_metrics(ds_test, lb, vref_dict; 
                               e_test = get_all_energies(ds_test),
                               f_test = get_all_forces(ds_test))
    coeff_norm = norm(lb.β)
    @show coeff_norm

    e_test_pred = get_all_energies_w_onebody(ds_test, lb, vref_dict)
    e_test_mae, e_test_rmse, e_test_rsq = calc_metrics(e_test, e_test_pred)
    @show e_test_mae, e_test_rmse


    f_test_pred = get_all_forces(ds_test,lb)
    f_test_mae, f_test_rmse, f_test_rsq = calc_metrics(f_test_pred, f_test)
    @show f_test_mae, f_test_rmse
    
    return e_test_mae, e_test_rmse, f_test_mae, f_test_rmse, coeff_norm
end

ds1_dict = load("../files/small_iso17_dataset.jld2")
ds_test  = ds1_dict["ds_test"]
println("length of ds test: $(length(ds_test))")

ds2_dict = load("../files/10K_train_configs.jld2")
ds_train = ds2_dict["ds_train"]

# computing average energy, to subtract from the dataset
e_train = get_all_energies(ds_train)
e_test  = get_all_energies(ds_test)
f_test  = get_all_forces(ds_test) #just computing this here, used later for error metrics

num_atoms = length(get_system(ds_train[1])) # all systems have same number of atoms, entire ISO17 database

avg_energy_per_atom = mean(vcat(e_train,e_test))/num_atoms
vref_dict = Dict(:H => avg_energy_per_atom,
                 :C => avg_energy_per_atom,
                 :O => avg_energy_per_atom)

adjust_energies(ds_train,vref_dict)


# (body order, poly degree)
pod_param_files = [#"4body_COH_param_set1.pod",
                   #"4body_COH_param_set2.pod",
                   #"4body_COH_param_set3.pod",
                   #"4body_COH_param_set4.pod",
                   #"4body_COH_param_set5.pod",
                   #"4body_COH_param_set6.pod",
                   #"4body_COH_param_set7.pod",
                   #"4body_COH_param_set8.pod",
                   #"4body_COH_param_set9.pod",
                   #"4body_COH_param_set10.pod",
                   #"4body_COH_param_set11.pod",
                   "4body_COH_param_set12.pod",
                   "4body_COH_param_set13.pod",
                   "4body_COH_param_set14.pod"]

ddict = Dict()
for pfile in pod_param_files
    println("param file : $(pfile)")

    lmp_pod = LAMMPS_POD(pfile, [:C,:O,:H])
    println("POD length: $(length(lmp_pod))")
   
#    lb = LBasisPotential(lmp_pod)
#    ###### fitting procedure
#    println("ooc fitting")
#    ws = [30.0, 1.0]
#    _AtWA, _AtWb = PotentialLearning.ooc_learn!(lb, ds_train;ws=ws,symmetrize=false, λ=0.01)
#
#    println("computing test descriptors")
#    edescr_test = compute_local_descriptors(ds_test, lmp_pod)
#    fdescr_test = compute_force_descriptors(ds_test, lmp_pod)
#    local_ds_test = DataSet(ds_test .+ edescr_test .+ fdescr_test)
#
#
#    e_test_mae, e_test_rmse, f_test_mae, f_test_rmse, coeff_norm  = compute_error_metrics(
#                                                                        local_ds_test,
#                                                                        lb,
#                                                                        vref_dict;
#                                                                        e_test = e_test,
#                                                                        f_test = f_test)
#    println("e_test_mae: $(e_test_mae), f_test_mae: $(f_test_mae)")
#    ddict[spec] = Dict("coeffs" => lb.β,
#                       "length" => length(ace),                        
#                       "norm"   => coeff_norm,
#                       "e_mae"  => e_test_mae,
#                       "e_rmse" => e_test_rmse,
#                       "f_mae"  => f_test_mae,
#                       "f_rmse" => f_test_rmse)
#
#    # I just need to be very cautious with memory usage                                                                                             
#    local_ds_test      = nothing
#    edescr_test  = nothing
#    fdescr_test  = nothing  
#    lb           = nothing 
#    _AtWA        = nothing
#    _AtWb        = nothing
#
#    GC.gc()
end

#save("summary_of_pod_basis_scan_iso17.jld2", Dict("ddict" => ddict))
