using Unitful
using UnitfulAtomic
using AtomsBase
using ElectronicStructure

a = 3.1964
c = 5.0511
bounding_box = [[a, 0.0, 0.0], 
                [-a*cosd(60), a*sind(60), 0.0], 
                [0.0, 0.0, c]]u"Å"

system = periodic_system([:Hf => [0, 0, 0],
                          :Hf => [1 / 3, 2 / 3, 1 / 2]], 
                         bounding_box, fractional=true)

qe_params = QeParameters(
    system=system,
    tstress=true,
    tprnfor=true,
    ecutwfc=40.0,
    conv_thr=1e-11,
    occupations="smearing",
    smearing="gaussian",
    degauss=0.01,
    mixing_beta=0.7,
    pseudopotentials=Dict("Hf"=>"Hf-sp.oncvpsp.upf"),
    kpts=(1, 1, 1),
    n_mpi_procs=2,
)

calc = QeCalculator()
state = QeState(qe_params)

calculate(calc, state)
e = auconvert(u"eV", energy(state))
println("Energy: $e")


