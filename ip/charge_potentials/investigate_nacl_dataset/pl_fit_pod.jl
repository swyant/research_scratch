using PotentialLearning, InteratomicPotentials 
using Unitful
using DelimitedFiles
using FileIO
using LinearAlgebra: norm

datapath =  "./train"
param_file = "../files/4body_NaCl_param.pod"
coeff_file = "../files/dummy_nacl_coefficients.pod"

# from subsampling_dpp.jl in PL.jl examples
function readext(path::String, ext::String)
    dir = readdir(path)
    substr = [split(f, ".") for f in dir]
    id = findall(x -> x[end] == ext, substr)
    return dir[id]
end

# from subsampling_dpp.jl in PL.jl examples
function concat_dataset(confs::Vector{DataSet})
    N = length(confs)
    confs_vec = [[confs[i][j] for j = 1:length(confs[i])] for i = 1:N]
    confs_all = reduce(vcat, confs_vec)
    return DataSet(confs_all)
end

file_arr = readext(datapath, "xyz")
confs_arr = [load_data(datapath*"/"*file, ExtXYZ(u"eV", u"Å")) for file in file_arr]

confs = concat_dataset(confs_arr)
ds = confs

lmp_pod = LAMMPS_POD(param_file, [:Na, :Cl])

# should make it so that I can create a dummy coeff file
# unfortunately need a way of figuring out how many basis functions there will be.
nacl_lbp = LBasisPotential(lmp_pod,coeff_file) 

ws, int = [1.0, 0.0], false
@time _AtWA, _AtWb  = PotentialLearning.ooc_learn!(nacl_lbp, ds; 
                                                   ws=ws,
                                                   reg_style=:scale_thresh,
                                                   #reg_style=:default,
                                                   symmetrize=true, 
                                                   eweight_normalized=nothing,
                                                   λ=1e-12)
                                                   #λ=0.01)

open("nacl_lbp_beta.txt", "w") do io
    writedlm(io, nacl_lbp.β)
end

natoms = [length(position(sys)) for sys in get_system.(ds)]
e_ref = get_all_energies(ds)
epa_ref = e_ref ./ natoms
f_ref = get_all_forces(ds)

@time begin
    e_descr = compute_local_descriptors(confs, lmp_pod)
    f_descr = compute_force_descriptors(confs, lmp_pod)
end
ds_test = DataSet(confs .+ e_descr .+ f_descr)

e_pred = get_all_energies(ds_test,nacl_lbp)
epa_pred = e_pred ./ natoms

@show e_mae, e_rmse, e_rsq = calc_metrics(epa_pred,epa_ref)

f_pred = get_all_forces(ds_test,nacl_lbp)

@show f_mae, f_rmse, f_rsq = calc_metrics(f_pred,f_ref)
