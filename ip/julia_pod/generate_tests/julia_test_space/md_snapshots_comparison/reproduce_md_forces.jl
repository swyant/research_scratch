using InteratomicPotentials, PotentialLearning
using Unitful
using Random
using LinearAlgebra: norm
using Statistics: mean 

#ref_xyz = "../files/Hf128_md-snapshots_default_4-body.xyz"
#param_file = "../files/sample_4body_hf_param.pod"
#coeff_file = "../files/ref_vsmall_lammps_compat_Hf_coefficients.pod"
#elems = [:Hf]

#ref_xyz = "../files/monoclinic_hfo2_md-snapshots_default_4-body.xyz"
#param_file = "../files/sample_4body_hfo2_param.pod"
#coeff_file = "../files/ref_vsmall_lammps_compat_HfO2_coefficients.pod"
#elems = [:Hf, :O]

#ref_xyz = "../files/small_Hf-snapshots_default_4-body.xyz"
#param_file = "../files/sample_4body_hf_param.pod"
#coeff_file = "../files/ref_vsmall_lammps_compat_Hf_coefficients.pod"
#elems = [:Hf]

ref_xyz = "../files/iso17_md-snapshots_default_4-body.xyz"
param_file = "../files/sample_4body_iso17_param.pod"
coeff_file = "../files/ref_1K_train_iso17_coefficients.pod"
elems = [:C, :O, :H]

configs = load_data(ref_xyz, ExtXYZ(u"eV", u"Å"))
ref_e = get_all_energies(configs)
ref_f = get_all_forces(configs)
lmp_pod = LAMMPS_POD(param_file, elems)

lb = LBasisPotential(lmp_pod, coeff_file)
edescr = compute_local_descriptors(configs,lmp_pod)
fdescr = compute_force_descriptors(configs,lmp_pod)
ds = DataSet(configs .+ edescr .+ fdescr)
julia_e = get_all_energies(ds,lb)
julia_f = get_all_forces(ds,lb)

@show maximum(abs.(julia_e .- ref_e))
@show maximum(abs.(julia_f .- ref_f))
#    lmp_pod = LAMMPS_POD(pfile, [:C,:O,:H])
#    println("POD length: $(length(lmp_pod))")
#   
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
#    ddict[pfile] = Dict("coeffs" => lb.β,
#                       "length" => length(lmp_pod),                        
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
#end
#
#save("summary_of_pod_basis_scan_iso17.jld2", Dict("ddict" => ddict))
