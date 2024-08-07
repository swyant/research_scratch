using AtomsIO

include("./v2_molly_interpot_utils.jl")

# load system with AtomsIO
sys = load_system(ExtxyzParser(), "dump_final.xyz")


# set up InteratomicPotentials LennardJones, based off of 10.1103/PhysRevB.54.340 
ϵ = 0.01032u"eV"
σ = 3.405u"Å"
rcut = 8.51u"Å"
species = [:Ar]
lj_p = InteratomicPotentials.LennardJones(ϵ,σ,rcut,species)

inter_lj = InteratomicPotentialInter(lj_p, u"eV", u"Å")
general_inters = (inter_lj,)

# Check force with regular sys obtained with AtomsIO
f_p = Molly.forces(inter_lj, sys)
#f_p = [uconvert.(u"eV/Å", fi) for fi in f_p]
#=
 [0.0019808268568929366 eV Å⁻¹, 0.020532665135586152 eV Å⁻¹, -0.01717748614964176 eV Å⁻¹]
 [0.0159924490724765 eV Å⁻¹, 0.0035814481937410225 eV Å⁻¹, -0.04837639800541863 eV Å⁻¹]
 [-0.02028532914817166 eV Å⁻¹, -0.0056662413044934225 eV Å⁻¹, 0.020174725791425038 eV Å⁻¹]
 [0.04269416032056045 eV Å⁻¹, 0.01709052715183096 eV Å⁻¹, -0.07347766990200989 eV Å⁻¹]
 [-0.05756693626161153 eV Å⁻¹, 0.07178852927421989 eV Å⁻¹, -0.0682536204422086 eV Å⁻¹]
 [0.027714826904755307 eV Å⁻¹, -0.03688926738500493 eV Å⁻¹, 0.06014495409852316 eV Å⁻¹]
 [0.022093629592778885 eV Å⁻¹, -0.0473875161499566 eV Å⁻¹, -0.07037312491157584 eV Å⁻¹]
 ⋮
 [-0.038654990750588154 eV Å⁻¹, 0.005408750407728867 eV Å⁻¹, -0.019720919488669258 eV Å⁻¹]
 [-0.0840473772875225 eV Å⁻¹, -0.013779025325309849 eV Å⁻¹, -0.10451457718812059 eV Å⁻¹]
 [-0.0036897502500228786 eV Å⁻¹, -0.03082684112393483 eV Å⁻¹, 0.02954487951879671 eV Å⁻¹]
 [0.050654375591296395 eV Å⁻¹, -0.037732010577144745 eV Å⁻¹, 0.10572291009547777 eV Å⁻¹]
 [-0.08844559990709276 eV Å⁻¹, -0.02123266565685698 eV Å⁻¹, -0.03853123003462404 eV Å⁻¹]
 [0.08315291873407642 eV Å⁻¹, -0.06355239283697076 eV Å⁻¹, -0.03609454189340165 eV Å⁻¹]
 [-0.05570071724763319 eV Å⁻¹, 0.02575751022503323 eV Å⁻¹, -0.040882285839913296 eV Å⁻¹]
 =#
mp = molly_params(sys)

m_sys = System(;mp...,
            general_inters=general_inters, 
            loggers=(force=ForceLogger(typeof(1.0u"eV/Å"), 1),),
            force_units=u"eV/Å",
            energy_units=u"eV",
            #loggers=(force=ForceLogger(Float32, 1),),
            #force_units=NoUnits,
            )

#### run zero simulation
simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

simulate!(m_sys,simulator,0)

m_sys.loggers.force.history[1]
#=
 [0.0019808268568929366 eV Å⁻¹, 0.020532665135586152 eV Å⁻¹, -0.01717748614964176 eV Å⁻¹]
 [0.0159924490724765 eV Å⁻¹, 0.0035814481937410225 eV Å⁻¹, -0.04837639800541863 eV Å⁻¹]
 [-0.02028532914817166 eV Å⁻¹, -0.0056662413044934225 eV Å⁻¹, 0.020174725791425038 eV Å⁻¹]
 [0.04269416032056045 eV Å⁻¹, 0.01709052715183096 eV Å⁻¹, -0.07347766990200989 eV Å⁻¹]
 [-0.05756693626161153 eV Å⁻¹, 0.07178852927421989 eV Å⁻¹, -0.0682536204422086 eV Å⁻¹]
 [0.027714826904755307 eV Å⁻¹, -0.03688926738500493 eV Å⁻¹, 0.06014495409852316 eV Å⁻¹]
 [0.022093629592778885 eV Å⁻¹, -0.0473875161499566 eV Å⁻¹, -0.07037312491157584 eV Å⁻¹]
 ⋮
 [-0.038654990750588154 eV Å⁻¹, 0.005408750407728867 eV Å⁻¹, -0.019720919488669258 eV Å⁻¹]
 [-0.0840473772875225 eV Å⁻¹, -0.013779025325309849 eV Å⁻¹, -0.10451457718812059 eV Å⁻¹]
 [-0.0036897502500228777 eV Å⁻¹, -0.03082684112393483 eV Å⁻¹, 0.02954487951879671 eV Å⁻¹]
 [0.050654375591296395 eV Å⁻¹, -0.037732010577144745 eV Å⁻¹, 0.10572291009547777 eV Å⁻¹]
 [-0.08844559990709276 eV Å⁻¹, -0.02123266565685698 eV Å⁻¹, -0.03853123003462404 eV Å⁻¹]
 [0.08315291873407642 eV Å⁻¹, -0.06355239283697076 eV Å⁻¹, -0.03609454189340165 eV Å⁻¹]
 [-0.05570071724763319 eV Å⁻¹, 0.02575751022503323 eV Å⁻¹, -0.040882285839913296 eV Å⁻¹]
=#
f_check = [f for f in m_sys.loggers.force.history[1]]

# Check the forces 
lines = readlines("./ref_forces")
f_ref = [parse.(Float64,split(li)) for li in lines]u"eV/Å"

fcomp_errs = []
for i in 1:length(f_ref)
    for j in 1:3
        err = f_check[i][j]-f_ref[i][j]
        push!(fcomp_errs,err)
    end
end

@show maximum(fcomp_errs)
#7.577272143066693e-15 eV Å⁻¹

sys = staticAtoms(sys)

m_sys2 = System(Molly.System(sys,
                             u"eV",
                             u"eV/Å"),
                general_inters=general_inters, 
                loggers=(force=ForceLogger(typeof(1.0u"eV/Å"), 1),),)

simulate!(m_sys2,simulator,0)
f_new = m_sys2.loggers.force.history[1]
#=
 [0.0019808268568929366 eV Å⁻¹, 0.020532665135586152 eV Å⁻¹, -0.01717748614964176 eV Å⁻¹]
 [0.0159924490724765 eV Å⁻¹, 0.0035814481937410225 eV Å⁻¹, -0.04837639800541863 eV Å⁻¹]
 [-0.02028532914817166 eV Å⁻¹, -0.0056662413044934225 eV Å⁻¹, 0.020174725791425038 eV Å⁻¹]
 [0.04269416032056045 eV Å⁻¹, 0.01709052715183096 eV Å⁻¹, -0.07347766990200989 eV Å⁻¹]
 [-0.05756693626161153 eV Å⁻¹, 0.07178852927421989 eV Å⁻¹, -0.0682536204422086 eV Å⁻¹]
 [0.027714826904755307 eV Å⁻¹, -0.03688926738500493 eV Å⁻¹, 0.06014495409852316 eV Å⁻¹]
 [0.022093629592778885 eV Å⁻¹, -0.0473875161499566 eV Å⁻¹, -0.07037312491157584 eV Å⁻¹]
 ⋮
 [-0.038654990750588154 eV Å⁻¹, 0.005408750407728867 eV Å⁻¹, -0.019720919488669258 eV Å⁻¹]
 [-0.0840473772875225 eV Å⁻¹, -0.013779025325309849 eV Å⁻¹, -0.10451457718812059 eV Å⁻¹]
 [-0.0036897502500228786 eV Å⁻¹, -0.03082684112393483 eV Å⁻¹, 0.02954487951879671 eV Å⁻¹]
 [0.050654375591296395 eV Å⁻¹, -0.037732010577144745 eV Å⁻¹, 0.10572291009547777 eV Å⁻¹]
 [-0.08844559990709276 eV Å⁻¹, -0.02123266565685698 eV Å⁻¹, -0.03853123003462404 eV Å⁻¹]
 [0.08315291873407642 eV Å⁻¹, -0.06355239283697076 eV Å⁻¹, -0.03609454189340165 eV Å⁻¹]
 [-0.05570071724763319 eV Å⁻¹, 0.02575751022503323 eV Å⁻¹, -0.040882285839913296 eV Å⁻¹]
=#