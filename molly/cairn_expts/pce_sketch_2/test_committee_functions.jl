using PotentialLearning, InteratomicPotentials 
using CSV, DataFrames
using AtomsIO
using AtomsBase
using Statistics

cd("./pce_sketch_2")
include("./committee_AL_utils.jl")

ace = ACE(  species            = [:Hf, :O],
            body_order         = 4, 
            polynomial_degree  = 10,
            wL                 = 1.5, 
            csp                = 1.0,
            r0                 = 2.15,
            rcutoff            = 5.0)

sys = load_system("./test_tetrag_hfo2.xyz")

# Define committee 
fit_folder = "./initial_fits/"
pots = [let 
           coeffs = vec(Matrix(CSV.read(fit_folder*"fit$(i)_coeffs.csv",DataFrame,header=false)))
           lbp = LBasisPotential(coeffs,[0.0],ace)
         end
         for i in 1:10]
my_cmte_pot = CommitteePotential(pots,1);

# Test things
potential_energy(sys,my_cmte_pot)
force(sys,my_cmte_pot)

compute_all_energies(sys,my_cmte_pot)
compute_all_forces(sys,my_cmte_pot)


my_cmte_pot2 = CommitteePotential(pots[1:3],1)
sys2 = load_system("./dummy_hfo2.xyz")

test_energy = CommitteeEnergy()

compute(test_energy,sys,my_cmte_pot2)