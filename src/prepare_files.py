"""
Runs setup.sh to prepare files for the sra_to_features_pipeline project.
"""

import configparser
import glob
import subprocess
import sys
from pathlib import Path

def main():
    """
    Run setup.sh to prepare files for the sra_to_features_pipeline project.
    Generates genome.sizes file based on the reference FASTA file.
    
    This script is intended to be run before the main pipeline execution.
    It ensures that all necessary files are in place and properly configured.
    """
    # Start by patching the CNVpytor genome file to fix the split error
    print("Patching CNVpytor genome file...")
    patch_cnvpytor_genome_file()

    # Start the setup process
    setup_script = Path("setup.sh")
    
    if not setup_script.exists():
        print(f"Error: {setup_script} does not exist.")
        sys.exit(1)
    
    print(f"Running {str(setup_script)} to prepare files...")
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


def patch_cnvpytor_genome_file():
    """
    Finds the cnvpytor/genome.py file in the virtual environment and 
    patches the 'split' error by converting a PosixPath object to a string.
    """
    # Use glob to find the genome.py file within the virtual environment
    file_path_pattern = f"{sys.prefix}/lib/python*/site-packages/cnvpytor/genome.py"
    genome_file = glob.glob(file_path_pattern)

    if not genome_file:
        print("Error: Could not find cnvpytor/genome.py. Patch failed.")
        sys.exit(1)

    genome_file_path = genome_file[0]
    print(f"Found cnvpytor/genome.py at: {genome_file_path}")

    # Read the file's content
    with open(genome_file_path, 'r') as f:
        lines = f.readlines()

    # Apply the patch to the specific line
    patched_lines = []
    found_patch_target = False
    for i, line in enumerate(lines):
        # Look for the line 'fn = res.split("/")[-1]' which is the source of the error.
        if "fn = res.split(\"/\")[-1]" in line:
            # Replace the line with the fixed version
            patched_lines.append(line.replace("res.split(\"/\")", "str(res).split(\"/\")"))
            found_patch_target = True
            print(f"Successfully patched line {i+1} in {genome_file_path}")
        else:
            patched_lines.append(line)

    if not found_patch_target:
        print("Error: The target line for patching was not found. The file may have already been patched or a different version is installed.")
        return None

    # Write the patched content back to the file
    with open(genome_file_path, 'w') as f:
        f.writelines(patched_lines)

    print("CNVpytor genome.py successfully patched.")

    return 0


if __name__ == "__main__":
    main()
    print("Setup complete. You can now run the sra_to_features_pipeline.")
