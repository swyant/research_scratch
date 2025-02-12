using AtomsBase # For now, 0.3.5, because it's what PL.jl is based off of. 
using Unitful, UnitfulAtomic
using StaticArrays
using InteratomicPotentials
using Printf

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

function get_all_energies(configs::Vector{<:AtomsBase.FlexibleSystem}; strip_units=true)
    if strip_units
        all_energies = [ustrip.(get_energy(config)) for config in configs]
    else
        all_energies = [get_energy(config) for config in configs]
    end
    return all_energies
end

function get_all_energies(configs::Vector{<:AtomsBase.FlexibleSystem}, lbp::LBasisPotential)
    all_energies = [potential_energy(config,lbp) for config in configs]
    return all_energies
end


function get_all_atomic_charges(configs::Vector{<:AtomsBase.FlexibleSystem}; strip_units=true)
    if strip_units
        all_atomic_charges = [ustrip.(get_atomic_charges(config)) for config in configs]
    else
        all_atomic_charges = [get_atomic_charges(config) for config in configs]
    end

    return reduce(vcat,all_atomic_charges)
end

# LBP should be LBasisChargeModel
function get_all_atomic_charges(configs::Vector{<:AtomsBase.FlexibleSystem}, lbp; strip_units=true)
    if strip_units
        all_atomic_charges = [atomic_charges(config,lbp,ustrip(get_total_charge(config)); with_units=false) for config in configs]
    else
        all_atomic_charges = [atomic_charges(config,lbp,ustrip(get_total_charge(config)); with_units=true) for config in configs]
    end

    return reduce(vcat,all_atomic_charges)
end

function calc_mae_rmse_mape(x_pred, x)
    x_mae = sum(abs.(x_pred .- x)) / length(x)
    x_rmse = sqrt(sum((x_pred .- x) .^ 2) / length(x))
    x_mape = sum(abs.( (x_pred .- x)./x )) / length(x)

    x_mae, x_rmse, x_mape
end


function write_extxyz(out_fname::String, configs::Vector{<:AtomsBase.FlexibleSystem})
    open(out_fname, "w") do outf
        for config in configs 
            num_atoms = length(config)
            # Write number of atoms
            println(outf, num_atoms)

            # Construct lattice string
            latvecs = ustrip.(bounding_box(config)) 
            lattice_str = @sprintf("Lattice=\"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\" ",
                                   latvecs[1]..., latvecs[2]..., latvecs[3]...)

            # Property string
            property_str = "Properties=species:S:1:pos:R:3:forces:R:3:atomic_charge:R:1 "

            # PBC and energy string
            pbc_energy_str = @sprintf("pbc=\"T T T\" energy=%.16f charge=%.16f", ustrip(get_energy(config)), ustrip(get_total_charge(config)))

            ## Combine info line
            #info_line = lattice_str * property_str * stress_str * pbc_energy_str
            info_line = lattice_str * property_str * pbc_energy_str
            println(outf, info_line)

            # Write atom data
            for atom in config.particles 
                @printf(outf, "%s %21.16f %21.16f %21.16f %21.16f %21.16f %21.16f %21.16f \n",
                        string(atomic_symbol(atom)),
                        ustrip.(position(atom))...,
                        ustrip.(atom.data[:force])...,
                        ustrip(atom.data[:atomic_charge]))
            end
        end
    end
end

function learn_only_energy(
    configs::Vector{<:AtomsBase.FlexibleSystem},
    basis::BasisSystem,
    ref_energies::Vector{Float64};
    λ::Float64 = 0.01,
    pbar::Bool = true,
)
    basis_size = length(basis)

    AtA = zeros(basis_size,basis_size)
    Atb = zeros(basis_size)

    if pbar
        iter = ProgressBar(configs)
    else
        iter = configs 
    end

    for (i,config) in enumerate(iter)
        ref_energy = ref_energies[i]

        global_descrs = reshape(sum(compute_local_descriptors(config,basis)),:,1)'

        A = [global_descrs;]
        b = [ref_energy;]

        AtA .+= A'*A
        Atb .+= A'*b
    end

    reg_matrix = λ*Diagonal(ones(size(AtA)[1]))
    AtA += reg_matrix

    β =  AtA \ Atb
    println("condition number of AtA: $(cond(AtA))")
    return β
end


# not sure why I'm getting zeros for all x components?
# make another auxiliary function takes in sys, outputs displaced sys

function make_displaced_structure(config::AtomsBase.AbstractSystem, i,alpha, ϵ; sys_dict=nothing)
    
    if isnothing(sys_dict)
        bbox      = bounding_box(config)
        bcs       = boundary_conditions(config)
        atom_syms = atomic_symbol(config)
        num_atoms = length(config)::Int64
    else
        bbox      = sys_dict["bbox"]
        bcs       = sys_dict["bcs"]
        atom_syms = sys_dict["atom_syms"]
        num_atoms = sys_dict["num_atoms"]
    end
   
    united_pos = position(config)
    length_unit = unit(eltype(united_pos[1])) # should be same unit for all positions
    disp_pos = ustrip.(united_pos)

    disp_pos[i] = setindex(disp_pos[i], disp_pos[i][alpha] + ϵ, alpha)
    disp_atoms = [AtomsBase.Atom(atom_syms[i], disp_pos[i]*length_unit) for i in 1:num_atoms]
    sys = AtomsBase.FlexibleSystem(disp_atoms,bbox, bcs) 

    sys
end

# -> 3N x N matrix of where each row is the derivative with respect to r_i_alpha
function finite_difference_charges(config::AtomsBase.FlexibleSystem, lbcm; ϵ = 0.001)
    
    total_charge = ustrip(get_total_charge(config))

    sys_dict = Dict{String,Any}()

    sys_dict["num_atoms"] = num_atoms = length(config)::Int64
    sys_dict["bbox"] = bounding_box(config)
    sys_dict["bcs"]  = boundary_conditions(config)
    sys_dict["atom_syms"] = atomic_symbol(config)
   
    
    charge_derivatives = Matrix{Float64}(undef, 3*num_atoms, num_atoms)

    for i in 1:num_atoms
        for alpha in 1:3
            forward_sys = make_displaced_structure(config,i,alpha, ϵ/2; sys_dict=sys_dict)
            forward_charges = atomic_charges(forward_sys, lbcm, total_charge; with_units=false)

            back_sys = make_displaced_structure(config,i,alpha, -ϵ/2; sys_dict=sys_dict)
            back_charges = atomic_charges(back_sys, lbcm, total_charge; with_units=false)

            charge_deriv = (forward_charges .- back_charges) ./ ϵ

            row_idx = 3*(i-1) + alpha
            charge_derivatives[row_idx,:]  .= charge_deriv
        end
   end

   charge_derivatives
end

# -> 3N x N matrix of where each row is the derivative with respect to r_i_alpha
function finite_difference_charges_old(config::AtomsBase.FlexibleSystem, lbcm; ϵ = 0.001)
    
    total_charge = ustrip(get_total_charge(config))

    num_atoms = length(config)::Int64
    bbox = bounding_box(config)
    bcs  = boundary_conditions(config)
    
    united_pos = position(config)
    length_unit = unit(eltype(united_pos[1])) # should be same unit for all positions

    base_pos = ustrip.(united_pos)
    atom_syms = atomic_symbol(config)
    
    charge_derivatives = Matrix{Float64}(undef, 3*num_atoms, num_atoms)

    for i in 1:num_atoms
        for alpha in 1:3
            disp_pos_forward = deepcopy(base_pos) # vector of SVector{3, Float64}
            disp_pos_forward[i] = setindex(disp_pos_forward[i], disp_pos_forward[i][alpha] + ϵ/2, alpha)
            forward_disp_atoms = [AtomsBase.Atom(atom_syms[i], disp_pos_forward[i]*length_unit) for i in 1:num_atoms]
            forward_sys = AtomsBase.FlexibleSystem(forward_disp_atoms,bbox, bcs) 
            forward_charges = atomic_charges(forward_sys, lbcm, total_charge; with_units=false)

            disp_pos_back = deepcopy(base_pos) # vector of SVector{3, Float64}
            disp_pos_back[i] = setindex(disp_pos_back[i], disp_pos_back[i][alpha] - ϵ/2, alpha)
            back_disp_atoms = [AtomsBase.Atom(atom_syms[i], disp_pos_back[i]*length_unit) for i in 1:num_atoms]
            back_sys = AtomsBase.FlexibleSystem(back_disp_atoms,bbox, bcs) 
            back_charges = atomic_charges(back_sys, lbcm, total_charge; with_units=false)

            #if i==1 && alpha==1
            #    println("disp_pos_forward:")
            #    display(disp_pos_forward)
            #    println("forward_charges:")
            #    display(forward_charges)
            #    println("disp_pos_back:")
            #    display(disp_pos_back)
            #    println("back_charges:")
            #    display(back_charges)
            #end

            charge_deriv = (forward_charges .- back_charges) ./ ϵ

            row_idx = 3*(i-1) + alpha
            charge_derivatives[row_idx,:]  .= charge_deriv
        end
   end

   charge_derivatives
end
