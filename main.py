#!/usr/bin/env python3


import os
import sys


if len(sys.argv) > 1:
    sra_file = sys.argv[1]
    print(f"The provided file with SRA IDs is: {sra_file}")
else:
    # Get the directory containing this script
    script_path = os.path.abspath(__file__)
    script_directory = os.path.dirname(script_path)
    # Define default SRA file name
    sra_file = os.path.join(script_directory, 'sra_table_selected.txt')
    print(f'No filename provided. Running with default file "{sra_file}".')

# Initialize dict-list of sras
sra_data = {}

l_sra = open(sra_file, 'r').read().split('\n')

for sra_id in l_sra:
    print(f'About to run pipeline.py with {sra_id}...')
    l = f'python3 pipeline.py {sra_id} /content/data'
    os.system(l)
