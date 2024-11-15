using PotentialLearning, InteratomicPotentials 
using Unitful
using CSV, DataFrames, Tables 

param_file = "../../files/sample_6body_hfo2_param.pod"
coeff_file = "../../files/dummy_coefficients.pod"

conf = load_data("./sample_rattled_HfO2.xyz", ExtXYZ(u"eV", u"Å"))
sys = get_system(conf[1])

lmp_pod = LAMMPS_POD(param_file, [:Hf,:O])
hfo2_lbp = LBasisPotential(lmp_pod,coeff_file)

fdescrs = compute_force_descriptors(sys,lmp_pod)
fdescrs_reshaped = [fdescrs[i][j] for i in eachindex(fdescrs) for j in 1:3]

fdescrs_reshaped = stack(fdescrs_reshaped; dims=1)

CSV.write("julia_computed_fdescrs.csv", DataFrame(Tables.table(fdescrs_reshaped)), header=false, delim=" "
