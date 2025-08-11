"""
Runs setup.sh to prepare files for the sra_to_features_pipeline project.
"""

import configparser
import sys
import subprocess
from pathlib import Path

def main():
    """
    Run setup.sh to prepare files for the sra_to_features_pipeline project.
    Generates genome.sizes file based on the reference FASTA file.
    
    This script is intended to be run before the main pipeline execution.
    It ensures that all necessary files are in place and properly configured.
    """
    setup_script = Path("setup.sh")
    
    if not setup_script.exists():
        print(f"Error: {setup_script} does not exist.")
        sys.exit(1)
    
    try:
        subprocess.run(["bash", str(setup_script)], check=True)
        print("Files prepared successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during setup: {e}")
        sys.exit(1)
    # Initialize configuration file parser
    config = configparser.ConfigParser()
    # Read the config.ini file
    config.read('config.ini')
    # Define base variables
    genome_sizes = config['Paths']['GENOME_SIZES']
    fasta_file = config['Paths']['REFERENCE_FASTA']

    # Creation of genome.sizes file
    print(f'Creating {genome_sizes} file...')
    with open(fasta_file, 'r') as f:
        dict_fasta = {}
        name = None
        length = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save records when finding a new one
                if name is not None:
                    dict_fasta[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
        # Save last record
        if name is not None:
            dict_fasta[name] = length
    # Save dict to file
    with open(genome_sizes, 'w') as f:
        for key, value in dict_fasta.items():
            str_out = f'{key}\t{value}\n'
            f.write(str_out)

    print(f'{genome_sizes} file created.')

    print('🎉 Setup complete.')


if __name__ == "__main__":
    main()
    print("Setup complete. You can now run the sra_to_features_pipeline.")
