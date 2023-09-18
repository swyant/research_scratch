using PotentialLearning

function subtract_peratom_e(config::Configuration, vref_dict)
    orig_e = get_energy(config)
    e_unit = orig_e.u
    new_e = get_values(orig_e)
    for atom in get_system(config).particles
        species = atom.atomic_symbol
        new_e -= vref_dict[species]
    end

    Energy(new_e,e_unit)
end