using Revise

using SteinMD
using Molly 
using PotentialLearning
using LinearAlgebra
using SpecialPolynomials

# plotting scripts - these are kept separate from the package for now
include("/Users/swyant/cesmix/dev/SteinMD.jl/src/makie/makie.jl")

ref = MullerBrown()
pce = PolynomialChaos(4, 2, ChebyshevU)

# define properties
atom_mass = 1.0u"g/mol"
boundary = RectangularBoundary(Inf*u"nm")
temp = 100.0u"K"

# define initial system
atoms = [Atom(mass=atom_mass)]
coords = [SVector(-0.8, 1.2)u"nm"] # initial position
sys = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    general_inters=(ref,),
    loggers=(coords=CoordinateLogger(100; dims=2),),
)

# define simulator
sim_langevin = OverdampedLangevin(
            dt=0.002u"ps",
            temperature=temp,
            friction=4.0u"ps^-1")

# run simulation - this will take a few seconds
simulate!(sys, sim_langevin, 5_000_000);

# subsample to obtain training data
coords_train = [sys.loggers.coords.history[i][1] for i = 2:2000:length(sys.loggers.coords.history)]
ntrain = length(coords_train)
atoms_train = [Atom(mass=atom_mass) for i in 1:ntrain]

sys_train = [System(
    atoms=[atoms_i],
    coords=[coords_i],
    boundary=boundary,
    general_inters=(ref,),
    # k = 1.0u"kJ * K^-1 * mol^-1",
) for (atoms_i, coords_i) in zip(atoms_train, coords_train)]


train_potential!(sys_train,ref,pce)


rbf = RBF(Euclidean(2), β=1.0, ℓ = 0.1)
# define simulator
sim_svgd = StochasticSVGD(
            dt=0.001u"ps",
            kernel=rbf,
            temperature=temp, #1.0u"K",
            friction=4.0u"ps^-1")
altrigger = TimeInterval(interval=200)

# define initial ensemble
ens0 = [System(
    atoms=[atoms_i],
    coords=[coords_i],
    boundary=boundary,
    general_inters=(pce,),
    # k = 1.0u"kJ * K^-1 * mol^-1",
    loggers=(
        coords=CoordinateLogger(1; dims=2),
        ksd=StepComponentLogger(1; dims=2),
        # trigger=TriggerLogger(altrigger, 1),
        params=TrainingLogger(),
    )
) for (atoms_i, coords_i) in zip(atoms_train, coords_train)]

ens = deepcopy(ens0) # I guess so that you constantly have access to the initialized ensemble
sys_final, alsteps = active_learn!(ens, sim_svgd, 2_000, sys_train, ref, altrigger)
# It's unstable