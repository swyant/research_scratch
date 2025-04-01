#! /usr/bin/env python
'''
using chatgpt w/ web search 3/24/25

Using the latest documentation for the atomic simulation environment package, write a python script that does the following:
1. reads in an xyz file with many different configurations, where the name of the file is an input argument to the script
2. randomly selects n configurations (where n is the second argument to the script) and writes this smaller set of configuration to an xyz file (the third argument of the script)

Critically, the input xyz file contains a system-level energy field and provides forces for each atoms, and so the output xyz should also have the appropriate energy and forces for each configuration
'''

import sys
import random
from ase.io import read, write

def main(input_file, num_configs, output_file):
    # Read all configurations from the input XYZ file
    configurations = read(input_file, index=':')
    
    # Ensure the number of configurations requested does not exceed available configurations
    if num_configs > len(configurations):
        raise ValueError(f"Requested number of configurations ({num_configs}) exceeds the total available ({len(configurations)}).")
    
    # Randomly select 'num_configs' configurations
    selected_configs = random.sample(configurations, num_configs)
    
    # Write the selected configurations to the output XYZ file
    write(output_file, selected_configs)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_xyz_file> <number_of_configurations> <output_xyz_file>")
        sys.exit(1)
    
    input_xyz_file = sys.argv[1]
    try:
        number_of_configurations = int(sys.argv[2])
    except ValueError:
        print("Error: The number of configurations must be an integer.")
        sys.exit(1)
    output_xyz_file = sys.argv[3]
    
    try:
        main(input_xyz_file, number_of_configurations, output_xyz_file)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
