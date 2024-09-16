using PotentialLearning, InteratomicPotentials 
using Unitful
using DelimitedFiles
using FileIO

datapath =  "../../files/small_lammps_compat_subset"
param_file = "../../files/sample_6body_hfo2_param.pod"
coeff_file = "../../files/dummy_coefficients.pod"

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

lmp_pod = LAMMPS_POD(param_file, [:Hf,:O])

@time e_descr = compute_local_descriptors(confs, lmp_pod)
@time f_descr = compute_force_descriptors(confs, lmp_pod)

ds = DataSet(confs .+ e_descr .+ f_descr)

save("small_lammps_compat_pod_descriptors.jld2", Dict("ds"=>ds))

# These coefficients are "dummy" coefficients, because they will be overwritten
# TODO: implement LBasisPotential constructor that instantiates dummy coeffs
#hfo2_lbp = LBasisPotential(lmp_pod,coeff_file)
#
#ws, int = [100.0, 2.0], false
#@time learn!(hfo2_lbp, ds, ws, int)
#
#open("hfo2_lbp_beta.txt", "w") do io
#    writedlm(io, hfo2_lbp.β)
#end
#
#natoms = [length(position(sys)) for sys in get_system.(ds)]
#e_ref = get_all_energies(ds)
#epa_ref = e_ref ./ natoms
#
#e_pred = get_all_energies(ds,hfo2_lbp)
#epa_pred = e_pred ./ natoms
#
#@show e_mae, e_rmse, e_rsq = calc_metrics(epa_pred,epa_ref)
