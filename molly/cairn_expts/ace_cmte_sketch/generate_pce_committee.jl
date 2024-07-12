using AtomsIO, Molly
using PotentialLearning, Cairn 
import StatsBase
using StaticArrays
using SpecialPolynomials 
using CSV, DataFrames 
using FileIO

include("./training_utils.jl")

ref = MullerBrownRot() 

sim_langevin =OverdampedLangevin(
                dt = 0.002u"ps", 
                temperature = 200.0u"K", 
                friction = 4.0u"ps^-1"
)

sys0 = System(ref, [0.5,0.5], loggers=(coords=CoordinateLogger(100;dims=2),))
simulate!(sys0, sim_langevin,1_000_000)

ids = StatsBase.sample(1:length(sys0.loggers.coords.history), 5000, replace=false)
full_trainset = [let 
                atom=[Molly.Atom(   mass=1.0u"g/mol", 
                                    σ=0.3u"nm", 
                                    ϵ=0.2u"kJ * mol^-1",)]
                coord = [SVector{2}(coord_container[1])] #Already have units!
                boundary = RectangularBoundary(Inf*u"nm")
                sys = Molly.System( atoms = atom,
                                    coords = coord, 
                                    boundary = boundary)
                end 
                for coord_container in sys0.loggers.coords.history[ids]] # for some reason, each entry in the history is a one-elem array of the SVector
save("full_trainset.jld2", Dict("full_trainset"=>full_trainset))

# PCE specification 
basisfam = Jacobi{0.5,0.5}
limits = [[-4.4,1.5],[-2,2]]
hq_order=20
lq_order=15

pce_hq = PolynomialChaos(hq_order,2, basisfam,xscl=limits)
pce_lq_template = PolynomialChaos(lq_order,2, basisfam,xscl=limits)

hq_lp = learn!(full_trainset, ref, pce_hq, [1000.,1],false;e_flag=true, f_flag=true)
CSV.write("./pce_fits/hq_init_coeffs.csv", DataFrame(Tables.table(hq_lp.β)),header=false)

pces = perform_and_record_pce_fits(full_trainset, pce_lq_template, ref, 10)