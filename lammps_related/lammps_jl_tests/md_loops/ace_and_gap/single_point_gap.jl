using Printf 

function extract_single_step_observables(lmp::LMP)
    atomids    = extract_atom(lmp, "id")
    raw_types  = extract_atom(lmp,"type")
    raw_pos    = extract_atom(lmp, "x")
    raw_forces = extract_atom(lmp, "f")
    raw_vels   = extract_atom(lmp, "v")

    pe = extract_compute(lmp,"pe",LAMMPS.API.LMP_STYLE_GLOBAL,LAMMPS.API.LMP_TYPE_SCALAR)
    
   
    sort_idxs = sortperm(atomids)
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


function get_gap_eandf(bbounds,atom_types,atom_pos,masses, tstep)
    num_types = size(masses)[1]
    num_atoms = size(atom_pos)[2]
    config_dict = LMP(["-screen", "none"]) do lmp
        command(lmp, "log none")

        command(lmp, "units          metal")
        command(lmp, "boundary       p p p")
        command(lmp, "atom_style     atomic")
        command(lmp, "neigh_modify   delay 0 every 1 check no") # probably not necessary for a single calc

        bbound_str = ""
        for i in 1:3
            for j in 1:2
                bbound_str = bbound_str * @sprintf " %.27f" bbounds[i][j]
            end
        end

        command(lmp, "region main block" * bbound_str)
        command(lmp, "create_box $(num_types) main ")
        
        for (type_idx,mass) in enumerate(masses)
            command(lmp, "mass $(type_idx) $(mass)")
        end

        # atom_types and atom_pos need to be in the same order, taken care of in extract function
        for j in 1:num_atoms
            xyz_str = ""
            for i in 1:3
                xyz_str = xyz_str * @sprintf " %.27f" atom_pos[i,j]
            end
            command(lmp, "create_atoms $(atom_types[j]) single" * xyz_str)
        end

        # obviously, this needs to parameterized rather than just hard-coded
        command(lmp, "pair_style  quip")
        command(lmp, """pair_coeff  * * gap.xml "Potential xml_label=GAP_2020_2_11_0_18_44_47_601" 72 8""")

        command(lmp, "dump           run_forces all custom 1 dump_gap.custom id type x y z fx fy fz")
        command(lmp, """dump_modify    run_forces append yes sort id format line "%4d %1d %32.27f %32.27f %32.27f %32.27f %32.27f %32.27f" """)

        command(lmp, "thermo       1")
        command(lmp, "thermo_style custom step temp pe ke etotal press")
        command(lmp,"compute pe all pe")
        
        command(lmp, "reset_timestep $(tstep)")
        command(lmp, "run 0")
        config_dict = extract_single_step_observables(lmp)
        return config_dict
    end
    return config_dict
end
