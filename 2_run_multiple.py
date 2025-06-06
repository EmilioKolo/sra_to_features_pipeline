#!/usr/bin/env python3

import os
import sys


venv_name = 'sra_to_vcf'

if len(sys.argv) > 1:
    sra_file = sys.argv[1]
    print(f"The provided file with SRA IDs is: {sra_file}")
else:
    sra_file = 'sra_table_selected.txt'
    print(f'No filename provided. Running with default file "{sra_file}".')

# Create a virtual environment
if not os.path.exists(venv_name):
    os.system(f'python3 -m venv ~/{venv_name}')
    print("Virtual environment created.")

# Activate the virtual environment
activate_script = f'~/{venv_name}/bin/activate'
os.system(f'source {activate_script}')

# Make install script executable
l = 'chmod +x ./0_install.sh'
os.system(l)

# Run install
l = './0_install.sh'
os.system(l)

# Initialize dict-list of sras
sra_data = {}

l_sra = open(sra_file, 'r').read().split('\n')

for sra_id in l_sra:
    print(f'About to run 1_main.py with {sra_id}...')
    l = f'python3 1_main.py {sra_id}'
    os.system(l)

