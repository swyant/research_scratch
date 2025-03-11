using PotentialLearning, InteratomicPotentials 
using Unitful
using InteratomicPotentials

param_file = "./sample_6body_hfo2_param.pod"
coeff_file = "./dummy_coefficients.pod"

lmp_pod = LAMMPS_POD(param_file, [:Hf,:O])

hfo2_lbp = LBasisPotential(lmp_pod,coeff_file)
