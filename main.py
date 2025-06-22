#!/usr/bin/env python3


import logging
import os
import sys


# Define logging print level
# Options: DEBUG, INFO, WARNING, ERROR or CRITICAL
logging.basicConfig(level=logging.INFO)

# Get system variables
if len(sys.argv) > 1:
    sra_file = sys.argv[1]
    sra_file = os.path.abspath(sra_file)
    logging.info(f"The provided file with SRA IDs is: {sra_file}")
    if len(sys.argv) > 2:
        output_dir = sys.argv[2].rstrip('/')
        # Use the absolute path for output_dir
        abs_output_dir = os.path.abspath(output_dir)
        logging.info(f'The output directory is "{abs_output_dir}"')
    else:
        output_dir = '../output'
        abs_output_dir = os.path.abspath(output_dir)
        w = f'Output directory not provided, using "{abs_output_dir}"'
        logging.info(w)
else:
    # Get the directory containing this script
    script_path = os.path.abspath(__file__)
    script_directory = os.path.dirname(script_path)
    # Define default SRA file name
    sra_file = os.path.join(script_directory, 'sra_table_selected.txt')
    w = f'No filename provided. Running with default file "{sra_file}".'
    logging.info(w)

# Define output folder inside docker and inside the computer
output_docker = '/content/data/output'
# Create the output directory if it does not exist
os.makedirs(abs_output_dir, exist_ok=True)

# Initialize dict-list of sras
sra_data = {}

l_sra = open(sra_file, 'r').read().split('\n')

for sra_id in l_sra:
    logging.info(f'About to run pipeline.py with {sra_id}...')
    l = 'docker run -u'
    l += f' -v {abs_output_dir}:{output_docker}'
    l += f' features-pipeline {sra_id} {output_docker}'
    os.system(l)
