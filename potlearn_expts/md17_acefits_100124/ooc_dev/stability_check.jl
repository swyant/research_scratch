using InteratomicPotentials, PotentialLearning
using Atomistic
using Unitful
using Molly 
using JLD2 
using Plots

#temperature=300u"K"
temperature=500u"K"

aspirin_xyzfile = "/Users/swyant/cesmix/datasets/revMD17/rmd17_aspirin.xyz"
full_ds = load_data(aspirin_xyzfile, ExtXYZ(u"eV", u"Å"))

raw_init_config = get_system(full_ds[1])

m_init_config = Molly.System(raw_init_config; 
                             energy_units=u"eV", 
                             force_units=u"eV/Å")

ace = ACE(species           = [:C,:O,:H],
      body_order        = 5,
      polynomial_degree = 9,
      wL                = 2.0,
      csp               = 1.0,
      r0                = 1.43,
      rcutoff           = 4.4 )

beta_dict = load("beta_sets.jld2")

lb = LBasisPotential(ace)
lb.β .= beta_dict["default_reg_β"]

ace_inter = InteratomicPotentialInter(lb,u"eV", u"Å")
general_inters = (ace_inter,)


m_sys = Molly.System(m_init_config;
                     general_inters=general_inters,
                     loggers=(energy=PotentialEnergyLogger(typeof(1.0u"eV"),1),)
                    ) 

simulator = NoseHoover( dt = 0.0001u"ps",
                        temperature = temperature,
                        remove_CM_motion = 1 )

random_velocities!(m_sys, temperature)

simulate!(m_sys, simulator, 5000)

#plot(m_sys.loggers.energy.history)
