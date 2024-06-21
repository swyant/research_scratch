using Cairn 
using Molly 
using SpecialPolynomials
#using LinearAlgebra 
#using PotentialLearning 

#include("./committee_AL_utils.jl")

ref = Himmelblau() 

temp = 100.0u"K"
sim_langevin = OverdampedLangevin(  dt=0.002u"ps",
                                    temperature=temp,
                                    friction=1.0u"ps^-1")
atom_mass = 1.0u"g/mol"
boundary = RectangularBoundary(Inf*u"nm")

init_sys = System(
           atoms=[Atom(mass=atom_mass)],
           coords=[SVector{2}([4.5,-2]) .* u"nm"],
           boundary=boundary,
           general_inters=(ref,),
           loggers=(
                  coords=CoordinateLogger(1000; dims=2),
                  energy=PotentialEnergyLogger(1000) #unused because training recomputes ref energy anyways
                  ),
           ) 

#simulate!(init_sys, sim_langevin, 1_000_000)
#
#coords_train = [init_sys.loggers.coords.history[i][1] for i in 1:length(init_sys.loggers.coords.history)]
#ntrain = length(coords_train)
#atoms_train = [Molly.Atom(mass=atom_mass) for i in 1:ntrain]
#
#full_train_set = [Molly.System(
#                  atoms=[atoms_i],
#                  coords=[coords_i],
#                  boundary=boundary,
#                  ) for (atoms_i, coords_i) in zip(atoms_train, coords_train)]
#
#
#limits = [[-6.25,6.25],[-5.75,5.75]]
#template_pce = PolynomialChaos(5, 2, Jacobi{0.5,0.5}, xscl=limits)
#cmte_pot, train_inds = train_committee!(template_pce,ref,full_train_set,10);
#
#primary_pce = deepcopy(template_pce)
#train_potential_e!(full_train_set,ref,primary_pce)
#
#trigger_lt = LowerThreshold(
