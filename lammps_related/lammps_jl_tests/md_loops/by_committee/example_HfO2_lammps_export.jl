using ACE1 
using ACE1pack
using ACEfit
using LinearAlgebra: I, Diagonal, pinv
using JuLIP
using CSV, DataFrames
using Random
using Interpolations
using OrderedCollections
using YAML

ds_path = "../../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz" 
raw_data = read_extxyz(ds_path)

train_idxs = vec(Matrix(CSV.read("../Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs   = vec(Matrix(CSV.read("../Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

rpib3 = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 3, 
            maxdeg  = 10,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 5.0,
            pin     = 0,
)

@show length(rpib3)

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))

vref = JuLIP.OneBody([:Hf => -2.70516846, :O => -0.01277342]...)
train_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data[train_idxs]]

val_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data[val_idxs]]

A3, Y3, W3 = ACEfit.assemble(train_data, rpib3)

solver_reg = ACEfit.QR(; lambda=1e-3)

results_reg3 = ACEfit.solve(solver_reg,A3,Y3)

@show results_reg3["C"]
CSV.write("N3_rcut5_maxdeg10_1e-3lambdaQR_fit_coeffs.csv", DataFrame(Tables.table(results_reg3["C"])),header=false)

pot_reg_tmp3 = JuLIP.MLIPs.combine(rpib3,results_reg3["C"])
pot_reg3 = JuLIP.MLIPs.SumIP(pot_reg_tmp3,vref)
ACE1pack.linear_errors(train_data,pot_reg3)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   5.961 │    0.268 │   0.000 │
│ crystal │   4.987 │    0.194 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.301 │    0.219 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│  amorph │   4.670 │    0.201 │   0.000 │
│ crystal │   3.900 │    0.142 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   4.133 │    0.160 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#

@show forces(pot_reg3,raw_data[395])
#=
108-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [0.2934028307365577, 2.2391463760462114, 3.580599125333118]
 [2.5019366314007585, 0.3079130993199866, -1.4811393668014048]
 [0.9859633875993222, 1.2623129729960136, -0.07248787421033143]
 [2.286952885143243, 1.810723932116241, 2.9778492198077755]
 [-0.5435944423851595, 4.861722008053475, -1.7159494358123766]
 [0.23237677472440918, 1.5272686727016698, 0.6949390029656103]
 [-0.8002883600232304, -2.2978908985803876, 0.37210545743083046]
 [-1.746547059810704, -0.14120515738752654, -0.5735341469154487]
 ⋮
 [0.1469309284411393, 0.14847947559621522, -2.137601068082883]
 [0.5138537348616117, 2.575780951551904, 2.7672776654489066]
 [0.6266532756672414, 4.539050978287776, -3.7780555685905277]
 [2.018868579233693, -0.695047855835776, -0.7356089824431966]
 [-2.4630114832334877, 1.2662416036869244, 2.4546481119625856]
 [-0.553912884394153, 0.8534391002040426, -1.436406495432046]
 [-0.20136050443277398, -1.7620995203815033, -0.7861252922056656]
 [4.826567066652775, 0.39094008889976917, -1.345187907057038]
=#

custom_export2lammps("Siviraman_GAP_HfO2_N3_rcut5_maxdeg10_regQR.yace", pot_reg3,rpib3)



######## UTILITY FUNCTIONS #########
function custom_export2lammps(fname, IP,rpib::RPIBasis)

    if length(IP.components) != 2
        throw("IP must have two components which are OneBody and ace")
    end

    ordered_components = []

    for target_type in [OneBody, PIPotential]
        did_not_find = true
        for i = 1:2
            if typeof(IP.components[i]) <: target_type
                push!(ordered_components, IP.components[i])
                did_not_find = false
            end
        end

        if did_not_find
            throw("IP must have two components which are OneBody and ace")
        end
    end

    V1 = ordered_components[1]
    V2 = ordered_components[2]

    species = collect(string.(chemical_symbol.(V2.pibasis.zlist.list)))
    species_dict = Dict(zip(collect(0:length(species)-1), species))
    reversed_species_dict = Dict(zip(species, collect(0:length(species)-1)))


    elements = Vector(undef, length(species))
    E0 = zeros(length(elements))

    for (index, element) in species_dict
        E0[index+1] = V1(Symbol(element))
        elements[index+1] = element
    end


    # Begin assembling data structure for YAML
    data = OrderedDict()
    data["elements"] = elements

    data["E0"] = E0

    # embeddings
    data["embeddings"] = Dict()
    for species_ind1 in sort(collect(keys(species_dict)))
        data["embeddings"][species_ind1] = Dict(
            "ndensity" => 1,
            "FS_parameters" => [1.0, 1.0],
            "npoti" => "FinnisSinclairShiftedScaled",
            "drho_core_cutoff" => 1.000000000000000000,
            "rho_core_cutoff" => 100000.000000000000000000)
    end

    # bonds
    data["bonds"] = OrderedDict()
    basis1p = deepcopy(rpib.pibasis.basis1p)
    radialsplines = ACE1.Splines.RadialSplines(basis1p.J; zlist = basis1p.zlist, nnodes = 10000)
    ranges, nodalvals, zlist = ACE1.Splines.export_splines(radialsplines)
    # compute spline derivatives
    # TODO: move this elsewhere
    nodalderivs = similar(nodalvals)
    for iz1 in 1:size(nodalvals,2), iz2 in 1:size(nodalvals,3)
        for i in 1:size(nodalvals,1)
            range = ranges[i,iz1,iz2]
            spl = radialsplines.splines[i,iz1,iz2]
            deriv(r) = Interpolations.gradient(spl,r)[1]
            nodalderivs[i,iz1,iz2] = deriv.(range)
        end
    end
    # ----- end section to move
    for iz1 in 1:size(nodalvals,2), iz2 in 1:size(nodalvals,3)
        data["bonds"][[iz1-1,iz2-1]] = OrderedDict{Any,Any}(
            "radbasename" => "ACE.jl",
            "rcut" => ranges[1,iz1,iz2][end],         # note hardcoded 1
            "nradial" => length(V2.pibasis.basis1p.J.J.A),
            "nbins" => length(ranges[1,iz1,iz2])-1)   # note hardcoded 1
        nodalvals_map = OrderedDict([i-1 => nodalvals[i,iz1,iz2] for i in 1:size(nodalvals,1)])
        data["bonds"][[iz1-1,iz2-1]]["splinenodalvals"] = nodalvals_map
        nodalderivs_map = OrderedDict([i-1 => nodalderivs[i,iz1,iz2] for i in 1:size(nodalvals,1)])
        data["bonds"][[iz1-1,iz2-1]]["splinenodalderivs"] = nodalderivs_map
    end

    functions, lmax = ACE1pack.export_ACE_functions(V2, species, reversed_species_dict)
    data["functions"] = functions
    data["lmax"] = lmax
    YAML.write_file(fname, data)
end