#!/usr/bin/env python
import json 
import re
import glob 
import math
from os.path import join
import random
import pprint
from lxml import etree as ET
import numpy as np 
import os 
import subprocess
import shutil

pp = pprint.PrettyPrinter(indent=4)

### Functions
def get_xml_files(xmldir, file_nums):
    filelist = []
    if file_nums  == '*':
        filelist = glob.glob(join(xmldir, 'vasprun*.xml'))

    elif isinstance(file_nums,list):
        for num_range in file_nums:
            single_num = re.match(r'^(\d+)$', num_range)
            number_range = re.match(r'^(\d+)\.\.(\d+)$', num_range)

            if single_num:
               filelist.append(join(xmldir,'vasprun{:d}.xml'.format(int(single_num[1]))))
            elif number_range and (int(number_range[1]) < int(number_range[2])):
               just_vaspruns = glob.glob(join(xmldir, 'vasprun*.xml'))
               filelist = filelist + [join(xmldir,f)  for f in just_vaspruns \
                                      if ( int(number_range[1]) <= int(re.findall(r'\d+', f)[-1]) <= int(number_range[2]) ) ]
            else:
                print("ERROR: incorrectly specifiec number range for xml files for " + xmldir)
    else: 
       print("ERROR: incorrectly specifiec file_nums dict value for " + xmldir) 

    filelist = sorted(filelist, key=lambda s: list(map(int, re.findall(r'\d+', s))))
    return filelist

def cull_xml_list(xmllist, cullfrac):
    cull_saved = math.floor(cullfrac*len(xmllist))
    shuffled_list = random.sample(xmllist, k=len(xmllist))

#    print("{:f} {:d} {:d}".format(cullfrac, len(xmllist),cull_saved))

    accepted_list = []
    rejected_list = []

    for i in range(cull_saved):
        accepted_list.append(shuffled_list[i])
    for j in range(cull_saved, len(shuffled_list)):
        rejected_list.append(shuffled_list[j])


    return accepted_list, rejected_list
    
def extract_xml_data(xml_file, species_list, first, **kwargs):
    if not first: 
        if "ref_dict" not in kwargs:
            print("need to pass reference dict")
            sys.exit(1)
        ref_dict = kwargs["ref_dict"]

    tree = ET.parse(xml_file)
    calc_dict = dict(zip(["latvecs","atom_types", "lammps_types", "crys_pos","cart_pos","forces","stress","energy"],[[],[],[],[],[],[],[],0.0]))

    for atom_species in tree.iterfind("./atominfo/array[@name='atoms']/set/rc/c[1]"):
        atom_species_text = atom_species.text.strip()
        try:
            relevant_species_dict = next(species_dict for species_dict in species_list if species_dict["atom_species"] == atom_species_text)  #This will have to change if I want to have two Al types
        except StopIteration:
            print("dictionary is not specified for atom species " +  atom_species + " in system " + xml_file) 
            sys.exit()

        species_deepmd_type = relevant_species_dict["deepmd_type"]
        species_lammps_type = relevant_species_dict["lammps_type"]
        calc_dict["atom_types"].append(species_deepmd_type)
        calc_dict["lammps_types"].append(species_lammps_type)
    

    # Check that it's the same "system"
    if not first:
      for i in range(len(calc_dict["atom_types"])):
          assert calc_dict["atom_types"][i] == ref_dict["atom_types"][i]
    
    calc = tree.find('calculation')

    for lv in calc.iterfind("./structure/crystal/varray[@name='basis']/v"):
        calc_dict["latvecs"].append([float(x) for x in re.split(r'\s+', lv.text.strip())])

    cry_to_cart = np.transpose(np.array(calc_dict["latvecs"]))

    for pos in calc.iterfind("./structure/varray[@name='positions']/v"):
        crys_pos = np.array([float(x) for x in re.split(r'\s+', pos.text.strip())])  #Need an additional check here if you want to have two Al types
        cart_pos = np.dot(crys_pos, cry_to_cart)
        calc_dict["crys_pos"].append(crys_pos)
        calc_dict["cart_pos"].append(cart_pos)
    
    if not first:
        assert ref_dict["num_atoms"] == len(calc_dict["crys_pos"])

    for force in calc.iterfind("./varray[@name='forces']/v"):
        calc_dict["forces"].append(np.array([float(x) for x in re.split(r'\s+', force.text.strip())]))

    temp_stress = []
    for stress in calc.iterfind("./varray[@name='stress']/v"):
        temp_stress.append([float(x) for x in re.split(r'\s+', stress.text.strip())])
    temp_stress = np.array(temp_stress)
    calc_dict["stress"] = temp_stress*1000  #kbar to bar


    calc_dict["energy"] = float(calc.find("./energy/i[@name='e_wo_entrp']").text.strip())

    calc_dict["xml_file"] = xml_file
    if first:
        first_ref_dict = {"atom_types" : calc_dict["atom_types"], "num_atoms" : len(calc_dict["crys_pos"])}
        return first_ref_dict, calc_dict
    else:
        return calc_dict


def gen_sys_rosetta(sys_name, sys_data, shuffled_full_train, shuffled_full_val, set_divisor):
    rosetta_list = []
#    for xml_data in [equil_info[1] for equil_info in sys_data["equil"]]:
    for equil_info in sys_data["equil"]:
        xml_data = equil_info[1]        
        xml_name = xml_data["xml_file"]

        if equil_info[0] == 0:
            rosetta_dict = { "xml_file" : xml_name, "isequil" : True, "fit_category" : "rejected", "sys_indices": [], "system_sets" : []} 
        elif equil_info[0] > 0:
            rosetta_dict = { "xml_file" : xml_name, "isequil" : True, "fit_category" : "", "sys_indices": [], "system_sets" : []} 
        else:
            print("Somethings wrong with an equil weight")
            sys.exit(1)

        rosetta_list.append(rosetta_dict)
    for xml_data in sys_data["regular"]:
        xml_name = xml_data["xml_file"]
        rosetta_dict = { "xml_file" : xml_name, "isequil" : False, "fit_category" : "", "sys_indices": [], "system_sets" : []}
        rosetta_list.append(rosetta_dict)
    for xml_data in sys_data["rejected"]:
        xml_name = xml_data["xml_file"]
        rosetta_dict = { "xml_file" : xml_name, "isequil" : False, "fit_category" : "rejected", "sys_indices": [], "system_sets" : []}
        rosetta_list.append(rosetta_dict)

    num_weighted_train = len(shuffled_full_train)
    

    fit_num = 1
    for xml_data in shuffled_full_train:
        xml_name = xml_data["xml_file"]
        
        #Fuck computational complexity, y'know what I'm sayin       
        rosetta_match = False
        for rosetta_dict in rosetta_list:
            if xml_name == rosetta_dict["xml_file"]:
                assert not rosetta_match
                rosetta_match = True
                
                if rosetta_dict["fit_category"] == "":
                    rosetta_dict["fit_category"] = "train"
                else:
                    assert rosetta_dict["fit_category"] == "train"

                rosetta_dict["sys_indices"].append(fit_num)

                set_num = fit_num // set_divisor
                data_location = "{:s}/set.{:03d}".format(sys_name, set_num)

                if data_location not in rosetta_dict["system_sets"]:
                    rosetta_dict["system_sets"].append(data_location)
        fit_num += 1
    
    assert (fit_num-1) == num_weighted_train

    for xml_data in shuffled_full_val:
        xml_name = xml_data["xml_file"]
        
        #Fuck computational complexity, y'know what I'm sayin       
        rosetta_match = False
        for rosetta_dict in rosetta_list:
            if xml_name == rosetta_dict["xml_file"]:
                assert not rosetta_match
                rosetta_match = True
                
                if rosetta_dict["fit_category"] == "":
                    rosetta_dict["fit_category"] = "val"
                else:
                    assert rosetta_dict["fit_category"] == "val"

                rosetta_dict["sys_indices"].append(fit_num)

                set_num = fit_num // set_divisor
                data_location = "{:s}/set.{:03d}".format(sys_name, set_num)

                if data_location not in rosetta_dict["system_sets"]:
                    rosetta_dict["system_sets"].append(data_location)
        fit_num += 1

    return rosetta_list
                
           


##Input
input_json = "sic-gan-aln_trainval_data.json"
num_train_sets = 4
write_virial = True

with open(input_json, "r") as inputjsonfile:
    input_data = json.load(inputjsonfile)

if not input_data["type"] == "train_val":
    print("input json is not train/val")
    sys.exit(1)

final_rosetta_list = []
for system in input_data["systems"]:
    system_data = {"equil"   : [], "regular" : [], "rejected" : []}  
    
    # Get data for EQUIL files and their weights
    first_xml = True
    for equil_dict in system["equil"]:
        if first_xml:
            ref_data, xml_data = extract_xml_data(equil_dict["xml_file"], input_data["species"], True) 
            first_xml = False
        else:
            xml_data = extract_xml_data(equil_dict["xml_file"], input_data["species"], False, ref_dict=ref_data)

        if "weight" in equil_dict:
            equilweight = equil_dict["weight"]
        else:
            equilweight = input_data["equil_weight"]

        system_data["equil"].append([equilweight, xml_data])
    
    # Determine files for regular data, cull if necessary, then extract data 
    for file_set in system["trainval"]:
        xml_list = get_xml_files(file_set["xml_dir"], file_set["file_nums"])
        cull_fraction = 1 
        if "cull" in file_set:
            #print("cull shows up for {:s}".format(file_set["xml_dir"]))
            cull_fraction = file_set["cull"]
        elif "cull_type" in file_set:
            #print("cull_type shows up for {:s}".format(file_set["xml_dir"]))
            culltype = file_set["cull_type"]
           
            if culltype not in input_data["cull_types"]:
                print("Uspecified cull type: {:s}".format(culltype))
                sys.exit(1)

            cull_fraction = input_data["cull_types"][culltype]

        if cull_fraction == 1:
            rejected = []
            accepted = xml_list
        elif 0 < cull_fraction < 1:
            accepted, rejected = cull_xml_list(xml_list, cull_fraction)
        elif cull_fraction == 0:
            rejected = xml_list
            accepted = []
        else:
            print("cull_fraction mis-specified for {:s}".format(file_set["xml_dir"]))
            sys.exit(1)

        for xmlfile in accepted:
            #Assuming there's always at least one equil file to build the reference dict for the system
            xml_data = extract_xml_data(xmlfile, input_data["species"], False, ref_dict=ref_data)
            system_data["regular"].append(xml_data)

        for xmlfile in rejected:
            xml_data = extract_xml_data(xmlfile, input_data["species"], False, ref_dict=ref_data)
            system_data["rejected"].append(xml_data)

    prelim_num_val = math.floor(input_data["val_frac"]*len(system_data["regular"])) #Excluding EQUIL files
    prelim_num_regular_train = len(system_data["regular"]) - prelim_num_val

    weighted_num_equil = 0
    for equil_data in system_data["equil"]:
        weighted_num_equil += equil_data[0]

    prelim_num_full_train = prelim_num_regular_train + weighted_num_equil 
    
    print("total regular: {:d},  prelim_val: {:d},  prelim_regular_train: {:d},  weighted_num_equil: {:d},  prelim_num_full_train: {:d}".format(len(system_data["regular"]), prelim_num_val, prelim_num_regular_train, weighted_num_equil, prelim_num_full_train))

    set_divisor = prelim_num_full_train // num_train_sets  #Want a set divisor that will split the train set into num_train_sets sets, with a remainder equal to num_val data
    bad_remainder = prelim_num_full_train % num_train_sets  #To fix num_val

    print("set_divisor : {:d},  bad_remainder: {:d}".format(set_divisor, bad_remainder))

    num_val = prelim_num_val + bad_remainder
    num_regular_train = prelim_num_regular_train - bad_remainder 
    num_full_train = num_regular_train + weighted_num_equil
    num_full_trainval = num_full_train + num_val

    print("final num_val: {:d}".format(num_val))
    assert num_full_trainval % set_divisor == num_val


    # Make the val set from the system_data["regular"], i.e. excluding EQUIL files
    shuffled_noequil_trainval = random.sample(system_data["regular"], k=len(system_data["regular"]))
    
    full_weighted_train_list = []
    for i in range(num_regular_train):
        full_weighted_train_list.append(shuffled_noequil_trainval[i])
    full_val_list = []
    for j in range(num_regular_train, len(shuffled_noequil_trainval)):  # = num_val
        full_val_list.append(shuffled_noequil_trainval[j])

    # randomly insert weighted EQUIL data into train list 
    for equil_dict in system_data["equil"]:
        for _ in range(equil_dict[0]):
            # https://stackoverflow.com/questions/2475518/python-how-to-append-elements-to-a-list-randomly
            full_weighted_train_list.insert(random.randrange(len(full_weighted_train_list)+1),equil_dict[1])


    assert len(full_weighted_train_list) == num_full_train
    assert len(full_val_list) == num_val
    assert ( (len(full_weighted_train_list) + len(full_val_list)) % set_divisor ) == num_val

    
    system_set_dir = os.path.join(".", format(system["name"]))
    os.mkdir(system_set_dir)

    # Make raw files

    box_rawfile = os.path.join(system_set_dir,"box.raw")
    with open(box_rawfile, "w") as boxraw:
        for train_datum in full_weighted_train_list:
            boxraw.write("{:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f}\n".format(*train_datum["latvecs"][0], *train_datum["latvecs"][1],*train_datum["latvecs"][2]))
        for val_datum in full_val_list:
            boxraw.write("{:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f}\n".format(*val_datum["latvecs"][0], *val_datum["latvecs"][1],*val_datum["latvecs"][2]))

    coord_rawfile = os.path.join(system_set_dir,"coord.raw")
    with open(coord_rawfile, "w") as coordraw:
        for train_datum in full_weighted_train_list:
            for cart_pos in train_datum["cart_pos"]:
                coordraw.write("{:13.8f} {:13.8f} {:13.8f} ".format(*cart_pos))
            coordraw.write("\n")
        for val_datum in full_val_list:
            for cart_pos in val_datum["cart_pos"]:
                coordraw.write("{:13.8f} {:13.8f} {:13.8f} ".format(*cart_pos))
            coordraw.write("\n")

    force_rawfile = os.path.join(system_set_dir,"force.raw")
    with open(force_rawfile, "w") as forceraw:
        for train_datum in full_weighted_train_list:
            for force in train_datum["forces"]:
                forceraw.write("{:13.8f} {:13.8f} {:13.8f} ".format(*force))
            forceraw.write("\n")
        for val_datum in full_val_list:
            for force in val_datum["forces"]:
                forceraw.write("{:13.8f} {:13.8f} {:13.8f} ".format(*force))
            forceraw.write("\n")

    energy_rawfile = os.path.join(system_set_dir,"energy.raw")
    with open(energy_rawfile, "w") as energyraw:
        for train_datum in full_weighted_train_list:
            energyraw.write("{:16.8f}\n".format(train_datum["energy"]))
        for val_datum in full_val_list:
            energyraw.write("{:16.8f}\n".format(val_datum["energy"]))

    if write_virial:
        virial_rawfile = os.path.join(system_set_dir,"virial.raw")
        with open(virial_rawfile, "w") as virialraw:
            for train_datum in full_weighted_train_list:
                virialraw.write("{:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f}\n".format(*train_datum["stress"][0], *train_datum["stress"][1],*train_datum["stress"][2]))
            for val_datum in full_val_list:
                virialraw.write("{:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f} {:13.8f}\n".format(*val_datum["stress"][0], *val_datum["stress"][1],*val_datum["stress"][2]))

    type_rawfile = os.path.join(system_set_dir, "type.raw")      
    with open(type_rawfile, "w") as typeraw:
        for atom_type_num in ref_data["atom_types"]:
            typeraw.write("{:d} ".format(atom_type_num))
        typeraw.write("\n")

    os.chdir(system_set_dir)
    rawtoset_process = subprocess.Popen(['../raw_to_set.sh',str(set_divisor)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    rawtoset_response, err = rawtoset_process.communicate()
    print(rawtoset_response.decode("utf-8").strip())
    
    orig_raw_folder = "./original_raw_files"
    os.mkdir(orig_raw_folder)
    shutil.move("./box.raw", orig_raw_folder)
    shutil.move("./energy.raw", orig_raw_folder)
    shutil.move("./force.raw", orig_raw_folder)
    shutil.move("./coord.raw", orig_raw_folder)
    if write_virial:
        shutil.move("./virial.raw", orig_raw_folder)
    os.chdir("../")

    # get system rosetta
    sys_rosetta = gen_sys_rosetta(system["name"], system_data, full_weighted_train_list, full_val_list, set_divisor)
    final_rosetta_list = final_rosetta_list + sys_rosetta

with open("rosetta.json", "w") as rosetta_json:
    json.dump(final_rosetta_list, rosetta_json)
