using AtomsIO 
using InteratomicPotentials

bad_sys = load_system("weird_blas_error/problem_O2.xyz")

lmp_pod = LAMMPS_POD("./sample_6body_hfo2_param.pod", [:Hf,:O])

compute_local_descriptors(bad_sys,lmp_pod)
