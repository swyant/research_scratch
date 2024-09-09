using InteratomicPotentials
using AtomsIO 
using LAMMPS 

hfo2_sys = load_system("../../files/sample_monoclinic_HfO2.xyz")
lmp_pod = LAMMPS_POD("../../files/sample_6body_hfo2_param.pod", [:Hf,:O])
hfo2_lbp = LBasisPotential(lmp_pod, "../../files/sample_6body_2elem_coeffs.pod")

InteratomicPotentials.setup_lammps_system!(hfo2_sys,lmp_pod)

atomtype_str = ""
for elem_symbol in lmp_pod.species_map
    atomtype_str = atomtype_str * " " * string(elem_symbol)
end
command(lmp_pod.lmp, """compute dd0  all podd/atom $(lmp_pod.param_file) "" "" $(atomtype_str)""")

command(lmp_pod.lmp, "run 0")

test_atomids = extract_atom(lmp_pod.lmp, "id"
test_types = extract_atom(lmp_pod.lmp, "type")
test_ld = extract_compute(lmp_pod.lmp,"ld", LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_ARRAY)

