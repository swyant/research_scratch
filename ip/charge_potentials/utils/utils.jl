import AtomsBase # For now, 0.3.5, because it's what PL.jl is based off of. 
using Unitful, UnitfulAtomic
using StaticArrays
using InteratomicPotentials

# right now, ignoring PL.jl framework, so no Configuration, no DataSet, no Energy, no Forces
# not well-sanitized
function read_behler_input(fname::String;
                           energy_units::Unitful.EnergyUnits = u"eV",
                           length_units::Unitful.LengthUnits=u"Å") # --> Vector{FlexibleSystem}, w/ data field populated w/ atomitc charges
    force_units = energy_units/length_units   

    atomic_systems = AtomsBase.FlexibleSystem[]
    bc = SVector{3,AtomsBase.BoundaryCondition}([AtomsBase.Periodic(), AtomsBase.Periodic(), AtomsBase.Periodic()])

    curr_box = SVector{3,Unitful.Length}[]
    curr_atoms = AtomsBase.Atom[]
    curr_ddict = Dict{Symbol, Any}(:total_energy   => 0.0*energy_units, 
                                   :total_charge   => 0.0u"e_au")
    in_sys = false
    count = 0 
    open(fname, "r") do f
        for line in eachline(f)
            tokens = split(strip(line))
            if tokens[1] == "begin"
                in_sys = true

            elseif tokens[1] == "end"
                if length(curr_box) == 0
                    curr_box = SVector{3}(SVector{3}([50.0,0.0,0.0]*length_units),
                                          SVector{3}([0.0,50.0,0.0]*length_units),
                                          SVector{3}([0.0,0.0,50.0]*length_units))
                elseif length(curr_box) == 3
                    curr_box = SVector{3}(curr_box)
                else
                    error("improper number of lattice vectors, should be 0 or 3")
                end

                push!(atomic_systems, AtomsBase.FlexibleSystem(curr_atoms, curr_box, bc, curr_ddict))

                curr_box = SVector{3,Unitful.Length}[]
                curr_atoms = AtomsBase.Atom[]
                curr_ddict = Dict{Symbol,Any}(:total_energy   => 0.0*energy_units, 
                                              :total_charge   => 0.0u"e_au")

                in_sys = false
                count += 1
            elseif !in_sys
                error("this probably shouldn't happen")

            elseif tokens[1] == "lattice"
                if length(tokens) != 4 # assuming 3D
                    error("Not enough components in lattice vector")
                end

                box_vec = parse.(Float64, tokens[2:4])*u"bohr" # units of the input file
                box_vec = uconvert.(length_units,box_vec)
                push!(curr_box, SVector{3}(box_vec)) 
            elseif tokens[1] == "atom"
                pos = uconvert.(length_units,parse.(Float64, tokens[2:4])*u"bohr")
                chem_species = Symbol(tokens[5])
                atomic_charge = parse(Float64, tokens[6])*u"e_au"
                forces = uconvert.(force_units,parse.(Float64, tokens[8:10])*u"hartree/bohr")
                
                atom= AtomsBase.Atom(chem_species, pos,
                                     atomic_charge = atomic_charge,
                                     force = forces)
                push!(curr_atoms,atom)
            elseif tokens[1] == "energy"
                sys_energy = uconvert(energy_units,parse(Float64,tokens[2])*u"hartree")
                curr_ddict[:total_energy] = sys_energy
            elseif tokens[1] == "charge"
                sys_charge = parse(Float64, tokens[2])*u"e_au"
                curr_ddict[:total_charge] = sys_charge
            else
                error("unknown keyword")
            end            
        end
    end
    atomic_systems
end


function get_atomic_charges(sys::AtomsBase.FlexibleSystem)
    atomic_charges = [particle[:atomic_charge] for particle in sys.particles]
    atomic_charges
end

function get_total_charge(sys::AtomsBase.FlexibleSystem)
    total_charge = sys[:total_charge]
    total_charge   
end

function get_energy(sys::AtomsBase.FlexibleSystem)
    total_energy = sys[:total_energy]
    total_energy
end

function get_forces(sys::AtomsBase.FlexibleSystem)
     forces = [particle[:force] for particle in sys.particles]
     forces  
end

function get_all_atomic_charges(configs::Vector{<:AtomsBase.FlexibleSystem}; strip_units=true)
    if strip_units
        all_atomic_charges = [ustrip.(get_atomic_charges(config)) for config in configs]
    else
        all_atomic_charges = [get_atomic_charges(config) for config in configs]
    end

    return reduce(vcat,all_atomic_charges)
end

function get_all_atomic_charges(configs::Vector{<:AtomsBase.FlexibleSystem}, lbp::LBasisPotential; strip_units=true)
    if strip_units
        all_atomic_charges = [atomic_charges(config,lbp; with_units=false) for config in configs]
    else
        all_atomic_charges = [atomic_charges(config,lbp; with_units=true) for config in configs]
    end

    return reduce(vcat,all_atomic_charges)
end

function calc_mae_rmse_mape(x_pred, x)
    x_mae = sum(abs.(x_pred .- x)) / length(x)
    x_rmse = sqrt(sum((x_pred .- x) .^ 2) / length(x))
    x_mape = sum(abs.( (x_pred .- x)./x )) / length(x)

    x_mae, x_rmse, x_mape
end
