using InteratomicPotentials
using AtomsIO

include("parse_pod_analysis_utils.jl")

data_folder= "../small_lammps_compat_subset"

lmp_pod=LAMMPS_POD("../sample_6body_hfo2_param.pod", [:Hf, :O])
hfo2_lbp = LBasisPotential(lmp_pod,"./ref_small_lammps_compat_HfO2_HIGHP_coefficients.pod")

header=[:config, :num_atoms, :volume, :energy, :DFT_energy, :energy_error, :force, :DFT_force, :force_error]
ddict = read_segmented_files("./ref_small_lammps_compat_HfO2_training_analysis.pod"; header=header, delim=' ', ignorerepeated=true, skipto=2)

#one_completed = false
for fname in keys(ddict)
#    if !one_completed
    @show xyz_fname =joinpath(data_folder, fname)
    config_set = load_trajectory(xyz_fname)
    my_energies = [ ]
    my_fqois = [ ]
    for config in config_set
        energy = potential_energy(config, hfo2_lbp)
        push!(my_energies, energy)

        pred_forces = force(config, hfo2_lbp)
        force_qoi = sqrt(sum(map(x-> x^2,reduce(hcat, pred_forces))))
        push!(my_fqois,force_qoi)
    end
    
    @show maximum(abs.(my_energies-ddict[fname][!,:energy]))
    @show maximum(abs.(my_fqois-ddict[fname][!,:force]))

#    global one_completed=true
#    end
end
