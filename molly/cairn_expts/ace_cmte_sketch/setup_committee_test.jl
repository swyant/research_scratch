using InteratomicPotentials 
using CSV, DataFrames
using AtomsIO
using Statistics
using Molly 
import LinearAlgebra: norm

#cd("./pce_sketch_2")
includet("./committee_AL_utils.jl")

ace = ACE(  species            = [:Hf, :O],
            body_order         = 4, 
            polynomial_degree  = 10,
            wL                 = 1.5, 
            csp                = 1.0,
            r0                 = 2.15,
            rcutoff            = 5.0)

raw_sys = load_system("./dummy_hfo2.xyz")

# Define committee 
fit_folder = "./initial_fits/"
pots = [let 
           coeffs = vec(Matrix(CSV.read(fit_folder*"fit$(i)_coeffs.csv",DataFrame,header=false)))
           lbp = LBasisPotential(coeffs,[0.0],ace)
         end
         for i in 1:10]
my_cmte_pot  = CommitteePotential(pots[1:3],1);
my_cmte_pot2 = CommitteePotential(pots[4:6],1);
my_cmte_pot3 = CommitteePotential(pots[7:10],1);

simulator = NoseHoover(dt=0.001u"ps",
                              temperature = 300u"K",
                              remove_CM_motion=1
                              )

m_sys = Molly.System(raw_sys, energy_units=u"eV", 
                              force_units=u"eV/Å")
m_sys = Molly.System(m_sys; 
                     general_inters = (my_cmte_pot,))                            

cmte_qoi1 = CommitteeEnergy(maximum)
trigger1  = CmteTrigger(cmte_qoi1,>,-13.0) # make sure to pass a float 
cmte_qoi2 = CommitteeFlatForces((cmte=std,coord_and_atom=mean))
trigger2  = CmteTrigger(cmte_qoi2,>,5.0; 
                       logger_spec=(:cmte_flat_forces,1))
cmte_qoi3 = CommitteeForces((coord=norm,cmte=std,atom=mean))
trigger3 = CmteTrigger(cmte_qoi3,>,8.0; 
                        logger_spec=(:avg_cmt_std_fmags, 1))

m_sys2 = Molly.System(m_sys;loggers=(coords=CoordinateLogger(10),
                                            force=ForceLogger(typeof(1.0u"eV/Å"),10),
                                            energy=PotentialEnergyLogger(typeof(1.0u"eV"),10))
                                            )

shared_trigger = SharedCmteTrigger(my_cmte_pot2,
                                   (trigger1, trigger2, trigger3))

shared_trigger2 = SharedCmteTrigger(my_cmte_pot2,
                                   (trigger1, trigger2, trigger3);
                                   energy_cache_field = :cmte_energies,
                                   force_cache_field = :cmte_forces)