#!/usr/bin/env python3

import os
import pandas as pd
import sys


if len(sys.argv) > 1:
    sra_file = sys.argv[1]
    print(f"The provided file with SRA IDs is: {sra_file}")
else:
    print("No filename provided.")
    raise ValueError


# Initialize dict-list of sras
sra_data = {}

l_sra = open(sra_file, 'r').read().split('\n')

for sra_id in l_sra:
    l = f'python3 1_main.py {sra_id}'
    os.system(l)

