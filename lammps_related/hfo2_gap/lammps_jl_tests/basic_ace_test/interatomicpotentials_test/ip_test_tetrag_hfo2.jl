using PotentialLearning
using InteratomicBasisPotentials 
using CSV, DataFrames
using Unitful, UnitfulAtomic
using AtomsBase
using StaticArrays

#### Custom functions (NEED TO MERGE THIS ASAP)
function fixed_load_data(file, extxyz::ExtXYZ; T=Float64)
    configs = Configuration[]
    open(file, "r") do io
        count = 1
        while !eof(io)
            # Read info line
            line = readline(io)
            num_atoms = parse(Int, line)

            line = readline(io)
            lattice_line = match(r"Lattice=\"(.*?)\" ", line).captures[1]
            lattice = parse.(Float64, split(lattice_line)) * extxyz.distance_units
            box = SVector{3}.([lattice[1:3], lattice[4:6], lattice[7:9]]) # MODIFIED: added SVector here
            energy = try
                energy_line = match(r"energy=(.*?) ", line).captures[1]
                energy = parse(Float64, energy_line)
                Energy(energy, extxyz.energy_units)
            catch
                Energy(NaN, extxyz.energy_units)
            end

            # try
            #     stress_line = match(r"stress=\"(.*?)\" ", line).captures[1]
            #     stress = parse.(Float64, split(stress_line))
            #     push!(stresses, SVector{6}(stress))
            # catch
            #     push!(stresses, SVector{6}( fill(NaN, (6,))))
            # end

            bc = []
            try
                bc_line = match(r"pbc=\"(.*?)\"", line).captures[1]
                bc = [t == "T" ? AtomsBase.Periodic() : DirichletZero() for t in split(bc_line)]
            catch
                bc = [DirichletZero(), DirichletZero(), DirichletZero()]
            end

            properties = match(r"Properties=(.*?) ", line).captures[1]
            properties = split(properties, ":")
            properties = [properties[i:i+2] for i = 1:3:(length(properties)-1)]
            atoms = Vector{AtomsBase.Atom}(undef, num_atoms)
            forces = Force[]
            for i = 1:num_atoms
                line = split(readline(io))
                line_count = 1
                position = 0.0
                element = 0.0
                data = Dict(())
                for prop in properties
                    if prop[1] == "species"
                        element = Symbol(line[line_count])
                        line_count += 1
                    elseif prop[1] == "pos"
                        position = SVector{3}(parse.(T, line[line_count:line_count+2]))
                        line_count += 3
                    elseif prop[1] == "move_mask"
                        ft = Symbol(line[line_count])
                        line_count += 1
                    elseif prop[1] == "tags"
                        ft = Symbol(line[line_count])
                        line_count += 1
                    elseif prop[1] == "forces" || prop[1] == "force" # MODIFIED: the property is force, not forces
                        push!(
                            forces,
                            Force(
                                parse.(T, line[line_count:line_count+2]),
                                extxyz.energy_units / extxyz.distance_units,
                            ),
                        )
                        line_count += 3
                    else
                        length = parse(Int, prop[3])
                        if length == 1
                            data = merge(data, Dict((Symbol(prop[1]) => line[line_count])))
                        else
                            data = merge(
                                data,
                                Dict((
                                    Symbol(prop[1]) =>
                                        line[line_count:line_count+length-1]
                                )),
                            )
                        end
                    end
                end
                atoms[i] = AtomsBase.Atom(element,position .* extxyz.distance_units) # MODIFIED: not including data
            end

            system = FlexibleSystem(atoms, box, bc)
            count += 1
            push!(configs, Configuration(system, energy, Forces(forces)))
        end
    end
    return DataSet(configs)
end

test_data_file = "./test_tetrag_hfo2.xyz"
test_ds = fixed_load_data(test_data_file, ExtXYZ(u"eV", u"Å"))

ace = ACE(  species            = [:Hf, :O],
            body_order         = 4, 
            polynomial_degree  = 10,
            wL                 = 1.5, 
            csp                = 1.0,
            r0                 = 2.15,
            rcutoff            = 5.0)

coeffs = vec(Matrix(CSV.read("./N3_rcut5_maxdeg10_1e-3lambdaQR_fit_coeffs.csv",DataFrame,header=false)))

lb = LBasisPotential(coeffs,ace)

f_descr_test = compute_force_descriptors(test_ds,ace)
test_ds = DataSet(test_ds .+ f_descr_test)

@show get_all_forces(test_ds,lb)
"""
0.29340283073895534
2.239146376044573
3.5805991253309912
2.5019366313966884
0.3079130993205581
-1.4811393667985158
0.9859633875991562
1.262312973001599
⋮
0.8534391002042412
-1.4364064954312425
-0.20136050443099407
-1.762099520374818
-0.7861252922086805
4.826567066645225
0.39094008889964016
-1.345187907054651
"""