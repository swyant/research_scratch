function extract_single_step_observables(lmp::LMP)
    atomids    = extract_atom(lmp, "id")
    raw_types  = extract_atom(lmp,"type")
    raw_pos    = extract_atom(lmp, "x")
    raw_forces = extract_atom(lmp, "f")
    raw_vels   = extract_atom(lmp, "v")

    pe = extract_compute(lmp,"pe",LAMMPS.API.LMP_STYLE_GLOBAL,LAMMPS.API.LMP_TYPE_SCALAR)
    
   
    sort_idxs = sortperm(atomids)
# Warning: If I didn't slice this, I'm not sure I'd have a valid memory reference afterwards. Should check!
    config_types  = raw_types[sort_idxs]
    config_pos    = raw_pos[:, sort_idxs] 
    config_forces = raw_forces[:, sort_idxs]
    config_vels   = raw_vels[:, sort_idxs]

    config_dict = Dict( "types"      => config_types,
                        "positions"  => config_pos,
                        "forces"     => config_forces,
                        "velocities" => config_vels,
                        "pe"         => pe)
    
    config_dict
end

function parse_pe_file(fname::String="./pe.dat"; config_list = Nothing)
    pe_data = CSV.read(fname, DataFrame, header=["tstep", "pe"], skipto=2, delim=" ")
    if config_list == Nothing 
        return pe_data
    else
        for (idx,pe_row) in enumerate(eachrow(pe_data))
            #A bit of a fragile assumption
            @assert config_list[idx]["timestep"] == pe_row.tstep
            config_list[idx]["pe"] = pe_row.pe
        end
        return config_list
    end
end
