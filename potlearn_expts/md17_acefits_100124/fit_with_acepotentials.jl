using CSV, DataFrames
using ACEpotentials

aspirin_dsfile  = "../../../../datasets/revMD17/rmd17_aspirin.xyz"
train_indexfile = "../../../../datasets/revMD17/splits/index_train_01.csv" 
test_indexfile  = "../../../../datasets/revMD17/splits/index_test_01.csv"

train_indices = vec(Matrix(CSV.read(train_indexfile,DataFrame,header=false)))
test_indices  = vec(Matrix(CSV.read(test_indexfile,DataFrame,header=false)))

raw_data = read_extxyz(aspirin_dsfile)

weights = Dict("default" => Dict("E" => 30.0, "F" => 1.0, "V" => 1.0))

# From https://github.com/ACEsuit/mace-jax/blob/4b899de2101c6e2085ee972aeac0e46a334fd9a0/configs/aspirin.gin#L38
Vref = OneBody(:H => -13.587222780835477,
               :C => -1029.4889999855063,
               :N => -1484.9814568572233,
               :O => -2041.9816003861047)

#data = [ AtomsData(at; energy_key = "energy", force_key = "forces", 
#                   weights=weights, v_ref=Vref) for at in raw_data]
#train_data = data[train_indices]
#test_data  = data[test_indices]

train_data = raw_data[train_indices]
test_data  = raw_data[test_indices]

#= acemodel
#
Confident:
rin  = 0.77 A
rcut = 4.4 A
(for aspirin), elements = [:C, :O, :H]
(correlation) order = 4

I think:
envelope default is (:x,2,2) which is pcut=pin=2
r0 should probably just be :bondlen

wL = 2.0
Dn = Dict("default" => 1.0)
Dl = Dict("default" => 2.0) #wL
Dd = Dict(1 => 20, 2 => 20, 3 => 20, 4 => 8)
DEG = ACE1.RPI.SparsePSHDegreeM(Dn, Dl, Dd)
D = DEG


Open Question:
:transform => (:agnesi, 2, 4), <~~ default is agnesi so I think that should be fine, 
I don't think there is a transform option that quite matches the paper, so I can play around with this a bit

pure2b = false


For Pair potential:
Confident:
pair_rcut = 5.5

I think:
- pair_rin is unknown, but it could either be rin or 0.0 (so pin=0),
actually doesn't look like it's used at all, I think there's just pair_rcut
pair_envelope =  (:r,2)

Open:
- no info about pair_degree so will use default 

will use defaults for everything

=#
elements = [:C,:O,:H]

Dn = Dict("default" => 1.0)
Dl = Dict("default" => 2.0) #wL
Dd = Dict(1 => 20, 2 => 20, 3 => 20, 4 => 8)
DEG = ACE1.RPI.SparsePSHDegreeM(Dn, Dl, Dd)

limited_Vref = OneBody(:H => -13.587222780835477,
                       :C => -1029.4889999855063,
                       :O => -2041.9816003861047)

model1 = acemodel(elements    = elements,
                  order       = 4,
                  rcut        = 4.4,
                  rin         = 0.77,
                  r0          = :bondlen,
                  totaldegree = Dict(1 => 20, 2 => 20, 3 => 20, 4 => 8),
                  wL          = 2.0,
                  transform   = (:agnesi,2,4), # default
                  envelope    = (:x, 2, 2),
                  pure2b      = false, # not default, can play around with this
                  pair_rcut   = 5.5,
                  Vref        = limited_Vref) # can there be more elements here than elements?
                  # defaults for everything else for the pair potential

solver =  ACEfit.LSQR(damp = 0.25, atol = 1e-6)

P = smoothness_prior(model1; p=1.5)

# acefit!() basically overwrites the existing ACEData, so have to pass in weights, keys, and vref all over
acefit!(model1, train_data; solver     = solver, 
                            prior      = P, 
                            weights    = weights,
                            energy_key = "energy", 
                            force_key  = "forces")
                            # will by default grab the model Vref
