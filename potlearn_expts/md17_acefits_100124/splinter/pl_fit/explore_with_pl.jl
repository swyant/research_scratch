using InteratomicPotentials, PotentialLearning
using CSV, DataFrames 
using Unitful
using LinearAlgebra 

include("./subtract_peratom_e.jl")

function adjust_energies(ds, vref_dict)
    for config in ds
        new_energy = subtract_peratom_e(config,vref_dict)
        config.data[Energy] = new_energy
    end
end

vref_dict = Dict(:H => -13.587222780835477,
                 :C => -1029.4889999855063,
                 :O => -2041.9816003861047)

ds_path = "../../../../../../datasets/revMD17/rmd17_aspirin.xyz"
ds = load_data(ds_path, ExtXYZ(u"eV", u"Å"))

train_indexfile = "../../../../../../datasets/revMD17/rmd17/splits/index_train_01.csv" 
test_indexfile  = "../../../../../../datasets/revMD17/rmd17/splits/index_test_01.csv"

train_idxs = vec(Matrix(CSV.read(train_indexfile,DataFrame,header=false)))
test_idxs  = vec(Matrix(CSV.read(test_indexfile,DataFrame,header=false)))

data_train  = ds[train_idxs]
data_test   = ds[test_idxs]


ace1 = ACE(species           = [:C,:O,:H],
          body_order        = 5,
          polynomial_degree = 9,
          wL                = 2.0,
          csp               = 1.0,
          r0                = 1.43,
          rcutoff           = 4.4 )

adjust_energies(data_train,vref_dict)
#adjust_energies(data_val,vref_dict)

println("computing training descriptors")
e_descr_train = compute_local_descriptors(data_train,ace1)
f_descr_train = compute_force_descriptors(data_train,ace1)
ds_train = DataSet(data_train .+ e_descr_train .+ f_descr_train)

#println("computing test descriptors")
#e_descr_test = compute_local_descriptors(data_test,ace1)
#f_descr_test = compute_force_descriptors(data_test,ace1)
#ds_test = DataSet(data_test .+ e_descr_test .+ f_descr_test)


lb = LBasisPotential(ace1)

#learn statement 
ws, int = [30/sqrt(21), 1.0], false

lp = PotentialLearning.LinearProblem(ds_train)

@views B_train = reduce(hcat, lp.B)'
@views dB_train = reduce(hcat, lp.dB)'
@views e_train = lp.e
@views f_train = reduce(vcat, lp.f)

@views A = [B_train; dB_train]
@views b = [e_train; f_train]

Q = Diagonal([ws[1] * ones(length(e_train));
              ws[2] * ones(length(f_train))])

βs = Vector{Float64}()
#learn!(lb, ds_train, ws, int)
