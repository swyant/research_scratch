using PotentialLearning, InteratomicPotentials 
using Unitful
using CSV, DataFrames, Tables 
using LAMMPS

param_file = "../../files/sample_6body_hfo2_param.pod"
coeff_file = "../../files/dummy_coefficients.pod"

conf = load_data("./sample_rattled_HfO2.xyz", ExtXYZ(u"eV", u"Å"))
sys = get_system(conf[1])

lmp_pod = LAMMPS_POD(param_file, [:Hf,:O])
hfo2_lbp = LBasisPotential(lmp_pod,coeff_file)

ref_fdescrs = Matrix(CSV.read("julia_computed_fdescrs.csv", DataFrame, header=false))

InteratomicPotentials.setup_lammps_system!(sys,lmp_pod)
LAMMPS.command(lmp_pod.lmp, "run 0")

#LAMMPS.API.lammps_extract_compute(lmp_pod.lmp, "gdd", LAMMPS.STYLE_GLOBAL, 2)
