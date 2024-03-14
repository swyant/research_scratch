using LAMMPS

lmp = LMP(["-screen", "none"])

command(lmp,"units          metal")
command(lmp,"boundary       p p p")
command(lmp,"atom_style     atomic")
command(lmp,"atom_modify    map yes")

command(lmp,"neighbor       0.5 bin")
command(lmp,"neigh_modify   every 1 delay 0 check yes")

#command(lmp,"read_data      test_monoclinic_hfo2_unit_DATA")
#command(lmp,"region main prism 0.0 5.1364755065543157 0.0 5.1934133375444702 0.0 5.2461598045212767 0.0 -0.88344296827430302 0.0")
command(lmp,"region main prism 0.0 2.0 0.0 2.0 0.0 2.0 0.0 0.0 0.0")
command(lmp,"create_box 2 main")

#command(lmp, "pair_style    zero 10.0")
#command(lmp, "pair_coeff    * *")

#command(lmp, """compute ld all pod/atom ../sample_6body_hfo2_param.pod "" ""  Hf O""")

command(lmp,"mass 1 178.486")
command(lmp,"mass 2 15.999")


command(lmp, "change_box all x final 0.000000000000000000000000000 5.136475506554315728635629057 y final 0.000000000000000000000000000 5.193413337544470209650171455 z final 0.000000000000000000000000000 5.246159804521276726063661044 xz final -0.883442968274303019882154331 boundary p p p")
command(lmp,"create_atoms 1 single 0.792968059999999974607476361 2.378919590000000194862650460 3.713793009999999839010342839")
command(lmp,"create_atoms 1 single 3.018340939999999861242940824 4.975626779999999804715571372 4.155442589999999825067789061")
command(lmp,"create_atoms 1 single 3.460063449999999818373908056 2.814494770000000034571030483 1.532364220000000054611177802")
command(lmp,"create_atoms 1 single 1.234690059999999922624169812 0.217787590000000003120916858 1.090714650000000007779021871")
command(lmp,"create_atoms 2 single 1.436610560000000091918082035 3.856244939999999843394107302 5.133774230000000216023181565")
command(lmp,"create_atoms 2 single 2.374698439999999965976940075 1.259537759999999950721871755 2.735461879999999901258433965")
command(lmp,"create_atoms 2 single 2.816420439999999913993633527 1.337169419999999941950363791 0.112382999999999996787458656")
command(lmp,"create_atoms 2 single 1.878332560000000039934775486 3.933876609999999995892494553 2.510695360000000153632981892")
command(lmp,"create_atoms 2 single 0.045238479999999997582804667 1.717767880000000024764972295 1.816496590000000077580466495")
command(lmp,"create_atoms 2 single 4.649515029999999882193151279 4.314475060000000361526417691 0.806582290000000035057325931")
command(lmp,"create_atoms 2 single 4.207792519999999925062184047 3.475646489999999921849394013 3.429660650000000199355554287")
command(lmp,"create_atoms 2 single -0.396484029999999987303738180 0.878939310000000029177158467 4.439574949999999908811787463")

command(lmp, "pair_style    zero 10.0")
command(lmp, "pair_coeff    * *")

command(lmp, """compute ld all pod/atom ../sample_6body_hfo2_param.pod "" ""  Hf O""")


#atomids = extract_atom(lmp, "id")
#sort_idxs = sortperm(atomids)
#pos = extract_atom(lmp,"x")
#pos = pos[:,sort_idxs]
#raw_types = extract_atom(lmp,"type")
#sorted_types = raw_types[sort_idxs,:]

#command(lmp,"pair_style     pod")
#command(lmp,"""pair_coeff     * * ../sample_6body_hfo2_param.pod ../sample_6body_2elem_coeffs.pod "" "" Hf O""")

command(lmp,"""compute dd all podd/atom ../sample_6body_hfo2_param.pod "" "" Hf O""")

### These commands have the same effect as the extract_compute command below, and cause the issue
#command(lmp,"variable sample_dd5 equal c_dd[5][49]")
#command(lmp,"""fix ddprint all print 1 "\${sample_dd5}" file sample_julia_dd.txt""")
#command(lmp, "thermo 1")
##command(lmp, "thermo_style custom step pe c_dd[5][49]")
#command(lmp, "thermo_style custom step pe v_sample_dd5")
#command(lmp, "dump 		     mydump_dd all custom 1 julia_dump_dd id type c_dd[*]")
#command(lmp, "dump_modify  mydump_dd sort id format float %20.15g")

#command(lmp,"rerun rerun.custom dump x y z")
command(lmp,"run 0")

#raw_ld = extract_compute(lmp,"ld", LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_ARRAY)'
#num_pod_types = 2 
#num_atoms = size(raw_ld)[1]
#num_perelem_ld = size(raw_ld)[2] + 1 # including 1-body terms
#total_num_ld = num_pod_types*(num_perelem_ld)
#final_dd = [[zeros(total_num_ld) for _ in 1:3] for __ in 1:num_atoms] 

### The extract_compute command is the issue
#raw_dd = extract_compute(lmp,"dd", LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_ARRAY)
#raw_dd = LAMMPS.extract_compute2(lmp,"dd", LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_ARRAY)
#sorted_dd = raw_dd[sort_idxs,:]


command(lmp, "delete_atoms group all")

command(lmp,"mass 1 178.486")
command(lmp,"mass 2 15.999")
command(lmp, "change_box all x final 0.000000000000000000000000000 5.136475506554315728635629057 y final 0.000000000000000000000000000 5.193413337544470209650171455 z final 0.000000000000000000000000000 5.246159804521276726063661044 xz final -0.883442968274303019882154331 boundary p p p")

command(lmp,"create_atoms 1 single 0.776203818919444787916006590 2.350166829158174763847455324 3.682797980824850103687140290")
command(lmp,"create_atoms 1 single 3.046909974681736965607115053 4.944150105556711416454618302 4.138358285486625653959436022")
command(lmp,"create_atoms 1 single 3.476827233049027299216504616 2.843247525261339170299379475 1.563359341695713888853447315")
command(lmp,"create_atoms 1 single 1.206121357063365095996232412 0.249264269623342554771383561 1.107798839749943642374319097")
command(lmp,"create_atoms 2 single 1.992388292737735655535402657 3.756099553796239387537525545 0.094680110732730332623319214")
command(lmp,"create_atoms 2 single 2.415659899710223434965428169 1.232091591716474976436757061 2.535088502166750856758881127")
command(lmp,"create_atoms 2 single 2.260642995344876116092791563 1.437314984791091765359283272 5.151477127466908001451884047")
command(lmp,"create_atoms 2 single 1.837371541104662497900790186 3.961322345404707068183824958 2.711068754608851882892395224")
command(lmp,"create_atoms 2 single -0.040098562735613671459411478 1.710173262166000407447086218 1.839590910655264810458220381")
command(lmp,"create_atoms 2 single 4.637992808257605403809975542 4.332359369278528937741157279 0.806657226745606470430516310")
command(lmp,"create_atoms 2 single 4.293129979308621280154056876 3.483241329365950100793725142 3.406566397691777936529433646")
command(lmp,"create_atoms 2 single -0.384961544867351790344400797 0.861055037948540880066161662 4.439500168075511510323849507")

#atomids = extract_atom(lmp, "id")
#sort_idxs = sortperm(atomids)
#pos = extract_atom(lmp,"x")
#pos = pos[:,sort_idxs]
#raw_types = extract_atom(lmp,"type")
#sorted_types = raw_types[sort_idxs,:]

command(lmp, "run 0")

#raw_ld = extract_compute(lmp,"ld", LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_ARRAY)'
#num_pod_types = 2 
#num_atoms = size(raw_ld)[1]
#num_perelem_ld = size(raw_ld)[2] + 1 # including 1-body terms
#total_num_ld = num_pod_types*(num_perelem_ld)
#final_dd = [[zeros(total_num_ld) for _ in 1:3] for __ in 1:num_atoms] 
raw_dd = extract_compute(lmp,"dd", LAMMPS.API.LMP_STYLE_ATOM,LAMMPS.API.LMP_TYPE_ARRAY)'
#sorted_dd = raw_dd[sort_idxs,:]

#@show sorted_dd[5,49]
