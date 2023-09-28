using ACE1, ACE1pack, ACEfit
using JuLIP
using CSV, DataFrames 

include("../../utilities/custom_export2lammps/custom_export2lammps.jl")

ds_path = "../../../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz"
raw_data = read_extxyz(ds_path)

rpib = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 2, 
            maxdeg  = 6,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 4.0,
            pin     = 0,
)

#rpib3 = ACE1.rpi_basis(;
#            species = [:Hf, :O],
#            N       = 3, 
#            maxdeg  = 10,
#            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
#            r0      = 2.15,
#            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
#            rcut    = 5.0,
#            pin     = 0,
#)

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))

vref = JuLIP.OneBody([:Hf => -2.70516846, :O => -0.01277342]...)
#train_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
#                        weights = weights, v_ref=vref) for at in raw_data[train_idxs]]
all_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data]

A3, Y3, W3 = ACEfit.assemble(all_data, rpib)
solver_reg = ACEfit.QR(; lambda=1e-3)
results_reg = ACEfit.solve(solver_reg,A3,Y3)

#@show results_reg3["C"]
#CSV.write("N3_rcut5_maxdeg10_1e-3lambdaQR_fit_coeffs.csv", DataFrame(Tables.table(results_reg3["C"])),header=false)

pot_reg_tmp = JuLIP.MLIPs.combine(rpib,results_reg["C"])
pot_reg = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
ACE1pack.linear_errors(all_data,pot_reg)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│ crystal │   0.006 │    0.286 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   0.006 │    0.286 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│ crystal │   0.006 │    0.230 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   0.006 │    0.230 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
custom_export2lammps("test.yace", pot_reg,rpib)

### Checking prior working ACE fit
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


prior_coeffs = vec(Matrix(CSV.read("N3_rcut5_maxdeg10_1e-3lambdaQR_fit_coeffs.csv",DataFrame,header=false)))
pot_reg_tmp3 = JuLIP.MLIPs.combine(rpib3,prior_coeffs)
pot_reg3 = JuLIP.MLIPs.SumIP(pot_reg_tmp3,vref)
ACE1pack.linear_errors(train_data,pot_reg3)
#=
[ Info: RMSE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│ crystal │   5.795 │    0.226 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.795 │    0.226 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
[ Info: MAE Table
┌─────────┬─────────┬──────────┬─────────┐
│    Type │ E [meV] │ F [eV/A] │ V [meV] │
├─────────┼─────────┼──────────┼─────────┤
│ crystal │   5.795 │    0.176 │   0.000 │
├─────────┼─────────┼──────────┼─────────┤
│     set │   5.795 │    0.176 │   0.000 │
└─────────┴─────────┴──────────┴─────────┘
=#
custom_export2lammps("test_prior.yace", pot_reg3,rpib3)

