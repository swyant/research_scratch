using PotentialLearning, InteratomicPotentials 
using Unitful
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

#inpath = "./quick_test_training/"
#inpath = "./lammps_compat_train_subset/"
inpath = "./small_lammps_compat_subset/"

file_arr = readext(inpath, "xyz")
confs_arr = [load_data(inpath*file, ExtXYZ(u"eV", u"Å")) for file in file_arr]
confs = concat_dataset(confs_arr)

lmp_pod = LAMMPS_POD("./sample_6body_hfo2_param.pod", [:Hf,:O])

#e_descr = compute_local_descriptors(confs, lmp_pod)
#f_descr = compute_force_descriptors(confs, lmp_pod)
#
#ds = DataSet(confs .+ e_descr .+ f_descr)
#
#hfo2_lbp1 = LBasisPotential(lmp_pod1,"./sample_6body_2elem_coeffs.pod")
#
#lb = LBasisPotential(basis)
#ws, int = [100.0, 2.0], false
#learn!(lb, ds_train, ws, int)
