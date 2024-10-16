using InteratomicPotentials, PotentialLearning
using Unitful
using Molly 

temperature=300u"K"

aspirin_xyzfile = "/Users/swyant/cesmix/datasets/revMD17/rmd17_aspirin.xyz"
full_ds = load_data(aspirin_xyzfile, ExtXYZ(u"eV", u"Å"))

init_config = get_system(full_ds[1])

ace = ACE(species           = [:C,:O,:H],
      body_order        = 5,
      polynomial_degree = 9,
      wL                = 2.0,
      csp               = 1.0,
      r0                = 1.43,
      rcutoff           = 4.4 )

lb = LBasisPotential(ace)

ace_inter = InteratomicPotentialInter(lb,u"eV", u"Å")
general_inters = (ace_inter,)


m_sys = Molly.System(init_config;
                     general_inters=general_inters,
                    ) 

simulator = NoseHoover( dt = 0.0001u"ps",
                        temperature = temperature,
                        remove_CM_motion = 1 )

random_velocities!(m_sys, temperature)

simulate!(m_sys, simulator, 1000)
