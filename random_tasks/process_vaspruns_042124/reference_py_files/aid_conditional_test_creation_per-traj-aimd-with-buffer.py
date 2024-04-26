#!/usr/bin/env python
import json 
import re
import glob 
import math
from os.path import join
import pprint
import os 
import sys

pp = pprint.PrettyPrinter(indent=4)

viable_sample_freqs = (25,75)

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

def get_traj_bounds(xmldir):
  ref_filename = os.path.join(xmldir, "reference_list")
  
  traj_bounds = {}  
  prior_traj_num = 0
  all_config_ids = []

  prior_vasprun_num = 0
  prior_step_num = 0

  #if the file doesn't exist, I'll just let the error be raised and the program end. 
  with open(ref_filename, "r") as ref_file:
    for line in ref_file:
      ref_line_tokens = re.split(r"\s+", line.strip())
      vaspxml_file =ref_line_tokens[0] 
      xml_origin = ref_line_tokens[-1]

      path_tokens = re.split(r"/", xml_origin)
      
      # check that aimd is a token, as every "test"-able should be aimd, and aimd should be a dedicated folder somewhere in the path
      aimd_exists = False
      aimd_token_idx = None
      for token in path_tokens:
        if re.match("aimd",token,flags=re.IGNORECASE):
          aimd_exists = True
          aimd_token_idx = path_tokens.index(token)
      assert aimd_exists

      traj_num = int(re.match(r"trajectory_(\d+)$", path_tokens[-2]).group(1))
      vasprun_num = int(re.match(r"vasprun(\d+)\.xml$",vaspxml_file).group(1))
      step_num = int(re.match("step_(\d+)$", path_tokens[-1]).group(1))
       
      assert vasprun_num - prior_vasprun_num == 1 # just a sanity check on the vasp ordering

      config_map = (vasprun_num,step_num)

      if traj_num > prior_traj_num:
        if len(traj_bounds) > 0:
          traj_bounds[prior_traj_num]["bounds"].append(prior_vasprun_num)
          traj_bounds[traj_num] = {"bounds"       : [vasprun_num],
                                   "config_ids" : [config_map]
                                  }
        else:
          traj_bounds[traj_num] = {"bounds"       : [vasprun_num],
                                   "config_ids" : [config_map]
                                  }

        prior_traj_num = traj_num
        assert step_num in viable_sample_freqs # i.e., first step of traj should correspond to original sample frequency (although not necessarily true, see my aimd_refinement script)

      elif traj_num == prior_traj_num:
        step_diff = step_num - prior_step_num
        assert step_diff in viable_sample_freqs  # want to make sure my assumption about a consistent sample frequency across trajectories is correct

        traj_bounds[traj_num]["config_ids"].append(config_map)
        
      else:
        print("no idea why traj num has gone down")
        sys.exit(1)

      all_config_ids.append(config_map)

      prior_vasprun_num = vasprun_num
      prior_step_num = step_num

    traj_bounds[traj_num]["bounds"].append(prior_vasprun_num) # i.e. final trajectory ends on last vasprun number

  return traj_bounds, all_config_ids


##Input
input_json = "prelim_sic-aln-gan_all_data_with-conditional.json"
buffer_window = 40
test_window = 40

#traj_modifications = { ("sic",200,1) : 9000,
#                       ("sic",200,2) : 9000,
#                       ("sic",300,1) : 9000,
#                       ("sic",300,2) : 9000,
#                       ("sic",400,1) : 9000,
#                       ("sic",400,2) : 9000,
#                       ("interface_sic-aln_aln-lp",200,1): 14000,
#                       ("interface_sic-aln_aln-lp",200,3): 14000,
#                      }

#traj_modifications = { ("interface_sic-aln_aln-lp",300,1): 10000,
#                       ("interface_sic-aln_aln-lp",200,3): "exclude",
#                       ("interface_sic-aln_aln-lp",300,4): 10000,
#                       ("interface_sic-aln_aln-lp",300,5): 10000,
#                      }

with open(input_json, "r") as inputjsonfile:
    input_data = json.load(inputjsonfile)

for system in input_data["systems"]:
    num_sys_data = 0
    print("SYSTEM: " + str(system["name"]))
    for file_set in system["trainval"]:
        xml_list = get_xml_files(file_set["xml_dir"], file_set["file_nums"])
        traj_bounds, all_config_ids = get_traj_bounds(file_set["xml_dir"])
        
        if "modify" in file_set:
          for modify_traj_key  in file_set["modify"]:
            bound_traj_key = int(modify_traj_key)

            if file_set["modify"][modify_traj_key] == "exclude":
              traj_bounds.pop(bound_traj_key)

            elif isinstance(file_set["modify"][modify_traj_key], int):
              cutoff_step = file_set["modify"][modify_traj_key]
              
              prev_vasprun_num = -1
              truncate_success = False
              for cnfg_map in traj_bounds[bound_traj_key]["config_ids"]:
                curr_step = cnfg_map[1]                
                if curr_step > cutoff_step:
                  assert not prev_vasprun_num == -1
                  traj_bounds[bound_traj_key]["bounds"][1] = prev_vasprun_num
                  truncate_success = True
                  break
                   
                prev_vasprun_num = cnfg_map[0]

              assert truncate_success

            else:
              print("invalid trajectory modifier")
              sys.exit(1)

        #for key in traj_bounds:
        #  print(str(key) +":" + str(traj_bounds[key]["bounds"]))

        if file_set["istest"]:

          trainval_bounds = []
          test_bounds = []
          for traj_key in traj_bounds:
            bounds = traj_bounds[traj_key]["bounds"]
            trainval_lower = bounds[0]
            test_upper = bounds[1]
            test_lower = test_upper - test_window + 1
            trainval_upper = test_upper - buffer_window -test_window
            
            # check that bounds match expected vasprun indices
            assert all_config_ids[trainval_lower-1][0] == trainval_lower
            assert all_config_ids[trainval_upper-1][0] == trainval_upper
            assert all_config_ids[test_lower-1][0] == test_lower
            assert all_config_ids[trainval_upper-1][0] == trainval_upper

            # check that the step difference is 1 ps for test window and buffer
            test_window_time_range = test_window*viable_sample_freqs[0]  # assuming that I'm only considering the first sample freq
            buffer_window_time_range = buffer_window*viable_sample_freqs[0]
            assert all_config_ids[test_upper-1][1] - all_config_ids[test_lower-2][1] ==  test_window_time_range  # want to use the config right before test_lower
            assert all_config_ids[test_upper-1][1] - all_config_ids[trainval_upper-1][1] == test_window_time_range + buffer_window_time_range

            trainval_bounds.append((trainval_lower,trainval_upper))
            test_bounds.append((test_lower,test_upper))
          
          assert len(trainval_bounds) == len(test_bounds)
          print(file_set["xml_dir"])
          bounds_string = "trainval: ["
          for trainval_bnd in trainval_bounds:
            bounds_string = bounds_string + "\"{:d}..{:d}\", ".format(*trainval_bnd)
          bounds_string = bounds_string + "]   | test: ["
          for test_bnd in test_bounds:
            bounds_string = bounds_string + "\"{:d}..{:d}\", ".format(*test_bnd)
          bounds_string = bounds_string + "]"
          print(bounds_string)

        else:
          print(file_set["xml_dir"])
          print("'*'  del")
