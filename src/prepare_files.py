"""
Runs setup.sh to prepare files for the sra_to_features_pipeline project.
"""

import sys
import subprocess
from pathlib import Path

def main():
    """
    Run setup.sh to prepare files for the sra_to_features_pipeline project.
    
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

if __name__ == "__main__":
    main()
    print("Setup complete. You can now run the sra_to_features_pipeline.")
