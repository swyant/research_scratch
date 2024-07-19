using Cairn, AtomisticQoIs
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

cmte_indices = [let 
                 indices = vec(Matrix(CSV.read(fit_folder*"lq_fit$(i)_train_idxs.csv",DataFrame,header=false)))
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

m_sys1 = Molly.System(sys0;
                     general_inters=(pce_hq,),
                     loggers=(coords=CoordinateLogger(10),
                           force=ForceLogger(10),
                           energy=PotentialEnergyLogger(10)))


trigger0 = CmteTrigger(cmte_qoi1,>,-13.0,my_cmte_pot3)
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

#baby_alroutine = ALRoutine(MullerBrownRot(), pce_hq, full_trainset,GreedySelector(),InefficientLearningProblem())
#baby_alroutine.trainset, new_sys = update_trainset!(baby_alroutine.ss,test_sys,baby_alroutine)
#baby_alroutine.mlip = retrain!(baby_alroutine.lp,test_sys,baby_alroutine);

greedy = GreedySelector()
ilp = InefficientLearningProblem()

#new_trigger0 = typeof(trigger0)(trigger0; cmte_pot=my_cmte_pot)
#new_sharedtrigger = typeof(shared_trigger2) 

#test_sys = deepcopy(m_sys);
#cmte_lp = SubsampleAppendCmteRetrain(InefficientLearningProblem(;ref=ref),cmte_indices);
#my_alroutine = ALRoutine(ref,pce_hq,full_trainset,(trigger0,),greedy,ilp,(cmte_lp,),Dict());
#initialize_al_cache!(my_alroutine)
#@show my_alroutine.cache
#my_alroutine.trainset, my_alroutine.cache[:trainset_changes] = update_trainset!(my_alroutine.ss,test_sys,my_alroutine);
#@show my_alroutine.cache[:trainset_changes]
#@show my_alroutine.triggers[1].cmte_pot.members[1].params[1:5]
#my_alroutine.triggers = update_triggers!(my_alroutine.triggers,my_alroutine.trigger_updates, test_sys, my_alroutine)
#@show my_alroutine.triggers[1].cmte_pot.members[1].params[1:5]

#cmte_lp2 = SubsampleAppendCmteRetrain(InefficientLearningProblem(;ref=ref), cmte_indices)
#my_alroutine2 = ALRoutine(ref, my_cmte_pot, full_trainset, (trigger1,), greedy, cmte_lp2, nothing, Dict()) # no trigger update because using sys.general_inters for the committee potential
#my_alroutine2.trainset, my_alroutine2.cache[:trainset_changes] = update_trainset!(my_alroutine2.ss,test_sys,my_alroutine2);
#@show my_alroutine2.mlip.members[1].params[1:5]
#my_alroutine2.mlip = retrain!(my_alroutine2.lp, test_sys, my_alroutine2)
#@show my_alroutine2.mlip.members[1].params[1:5]


# use grid to define uniform quadrature points
coords_eval = potential_grid_2d(ref,limits,0.04,cutoff=800)
sys_eval = Ensemble(ref,coords_eval)
ζ = [ustrip.(Vector(coords)) for coords in coords_eval]
GQint = GaussQuadrature(ζ,ones(length(ζ)) ./length(ζ))
error_spec = Simple2DPotErrors(GQint,true)

aldata_spec = DefaultALDataSpec(error_spec,true,true,true,true)

#test_sys = deepcopy(m_sys1);
#cmte_lp = SubsampleAppendCmteRetrain(InefficientLearningProblem(;ref=ref),cmte_indices);
#my_alroutine = ALRoutine(ref,pce_hq,full_trainset,(shared_trigger2,),greedy,ilp,(cmte_lp,),aldata_spec,Dict());
#initialize_al_cache!(my_alroutine)
#test_sys = initialize_triggers(my_alroutine.triggers,test_sys)
#check log
#aldata = initialize_al_record(my_alroutine.aldata_spec)
#check aldata
#step_n = 10
#trigger_activated(my_alroutine.triggers,test_sys,my_alroutine, step_n)
#check logs
#my_alroutine.cache[:trigger_step] = step_n
#al.trainset, al.cache[:trainset_changes] = update_trainset!(my_alroutine.ss,test_sys,my_alroutine);
#@show my_alroutine.mlip.params
#my_alroutine.mlip = retrain!(my_alroutine.lp, test_sys, my_alroutine)
#@show my_alroutine.mlip.params
#@show my_alroutine.triggers[1].cmte_pot.members[1].params[1:5]
#@show my_alroutine.trigger_updates[1].cmte_indices[1][end-10:end]
#my_alroutine.triggers = update_triggers!(my_alroutine.triggers, my_alroutine.trigger_updates, test_sys, my_alroutine)
#@show my_alroutine.triggers[1].cmte_pot.members[1].params[1:5]
#@show my_alroutine.trigger_updates[1].cmte_indices[1][end-10:end]
#@show aldata 
#record_al_record!(my_alroutine.aldata_spec, aldata, test_sys, my_alroutine)
#@show aldata
#perstep_reset!(my_alroutine.triggers,test_sys)

