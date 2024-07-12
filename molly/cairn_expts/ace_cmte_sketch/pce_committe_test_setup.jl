using Cairn 
using SpecialPolynomials
using CSV, DataFrames
using Molly
using Statistics 
import LinearAlgebra: norm
using FileIO

includet("./committee_AL_utils.jl")

### reference potential and initial system
ref = MullerBrownRot()
sys0 = System(ref, [0.5,0.5])

### PCE specification 
basisfam = Jacobi{0.5,0.5}
limits = [[-4.4,1.5],[-2,2]]
hq_order=20
lq_order=15

### construct relevant pces
pce_hq = PolynomialChaos(hq_order,2, basisfam,xscl=limits)
pce_lq_template = PolynomialChaos(lq_order,2, basisfam,xscl=limits)

fit_folder ="./pce_fits/"
hq_coeffs = vec(Matrix(CSV.read(fit_folder*"hq_init_coeffs.csv",DataFrame,header=false)))
pce_hq.params = hq_coeffs

pots =  [let 
            coeffs = vec(Matrix(CSV.read(fit_folder*"lq_fit$(i)_coeffs.csv",DataFrame,header=false)))
            pce = deepcopy(pce_lq_template)
            pce.params = coeffs
            pce
        end
        for i in 1:10]

### setup committee potential
my_cmte_pot  = CommitteePotential(pots[1:3],1; 
                                  energy_units=pce_lq_template.energy_units,
                                  force_units=pce_lq_template.force_units)
my_cmte_pot2 = CommitteePotential(pots[4:6],1; 
                                  energy_units=pce_lq_template.energy_units,
                                  force_units=pce_lq_template.force_units)
my_cmte_pot3 = CommitteePotential(pots[6:10],1; 
                                  energy_units=pce_lq_template.energy_units,
                                  force_units=pce_lq_template.force_units)

### for simulate!() tests 
sim = OverdampedLangevin(dt=0.002u"ps",
                         temperature = 300.0u"K",
                        friction=4.0u"ps^-1")
m_sys0= Molly.System(sys0;
                     general_inters=(my_cmte_pot,),
                     loggers=(energy=PotentialEnergyLogger(1),))

#simulate!(m_sys0 sim, 10)

### setting up cmte_qois 
cmte_qoi1 = CommitteeEnergy(maximum)
#compute(cmte_qoi1,m_sys0,my_cmte_pot)
cmte_qoi2 = CommitteeFlatForces((cmte=std,coord_and_atom=mean))
#compute(cmte_qoi2,m_sys0,my_cmte_pot)
cmte_qoi3 = CommitteeForces((coord=norm, cmte=std, atom=mean))
#compute(cmte_qoi3,m_sys0,my_cmte_pot)


### Setting up trigger tests
m_sys = Molly.System(sys0;
                     general_inters=(my_cmte_pot,),
                     loggers=(coords=CoordinateLogger(10),
                           force=ForceLogger(10),
                           energy=PotentialEnergyLogger(10)))

trigger1 = CmteTrigger(cmte_qoi1,>,-13.0)
trigger2 = CmteTrigger(cmte_qoi2,>,5.0;
                       logger_spec=(:cmte_flat_forces,1))
trigger3 = CmteTrigger(cmte_qoi3,>,8.0;
                        logger_spec=(:avg_cmte_std_fmags,1))

shared_trigger = SharedCmteTrigger(my_cmte_pot2,
(trigger1, trigger2, trigger3))

shared_trigger2 = SharedCmteTrigger(my_cmte_pot2,
                        (trigger1, trigger2, trigger3);
                        energy_cache_field = :cmte_energies,
                        force_cache_field = :cmte_forces)

### Some Basic CmteTrigger and SharedCmteTrigger Testing 

#test_sys = deepcopy(m_sys)
#test_sys = initialize_triggers((trigger1,),test_sys)
#trigger_activated(trigger1,test_sys)
#perstep_reset!((trigger1,)test_sys)
#
#test_sys = deepcopy(m_sys)
#test_sys = initialize_triggers((trigger2,),test_sys)
#trigger_activated(trigger2,test_sys)
#perstep_reset!((trigger2,)test_sys)
#
#test_sys = deepcopy(m_sys)
#test_sys = initialize_triggers((shared_trigger,),test_sys)
#trigger_activated!(shared_trigger,test_sys)
#perstep_reset!((shared_trigger,),test_sys)
#
#test_sys = deepcopy(m_sys)
#test_sys = initialize_triggers((shared_trigger2,),test_sys)
#trigger_activated!(shared_trigger2,test_sys)
#perstep_reset!((shared_trigger2,),test_sys)


### loading in the full_trainset 

fdict = load("full_trainset.jld2")
full_trainset = fdict["full_trainset"]

#baby_alroutine = ALRoutine(full_trainset)
#test_sys = deepcopy(m_sys)
#@show position(test_sys)
#baby_alroutine.trainset, new_sys = update_trainset!(GreedySelector(),test_sys,baby_alroutine)
#@test position(new_sys) == [0.5,0.5]
#new_stripped_positions = [ustrip.(position(config)[1]) for config in baby_alroutine.trainset]
#@test any(map(x->x≈[0.5,0.5],new_stripped_positions)) == true


#baby_alroutine = ALRoutine(MullerBrownRot(), pce_hq, full_trainset)
#baby_alroutine.trainset, new_sys = update_trainset!(GreedySelector(),test_sys,baby_alroutine)
#@show old_params = baby_alroutine.mlip.params
#baby_alroutine.mlip = retrain!(InefficientLearningProblem(),test_sys,baby_alroutine)
#@show baby_alroutine.mlip.params
