#!/bin/bash

# ==============================================================================
# BASH TEST SCRIPT FOR SETUP.SH
#
# To run this script:
# 1. Ensure this file is in the same directory as your setup.sh script.
# 2. Make sure both scripts are executable:
#    chmod +x setup.sh test_setup.sh
# 3. Run the test script:
#    ./test_setup.sh
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# Define a temporary directory for the test.
TEMP_DIR=$(mktemp -d -t test-setup-XXXXXX)
echo "Using temporary directory: $TEMP_DIR"

# Define a function to clean up the temporary directory
cleanup() {
    echo "Cleaning up temporary directory..."
    rm -rf "$TEMP_DIR"
    echo "Cleanup complete."
}

# Register the cleanup function to be called on EXIT
trap cleanup EXIT

# ==============================================================================
# SETUP FUNCTIONS
# ==============================================================================

# Function to create mock files and directories needed by setup.sh
setup_mocks() {
    echo "Setting up mock environment..."

    # Create a mock config.ini file with test values
    cat <<EOF > "$TEMP_DIR/config.ini"
[Paths]
BASE_DIR = $TEMP_DIR/test_pipeline_base
KRAKEN_DB = $TEMP_DIR/test_pipeline_base/install/kraken_db
SNPEFF_DIR = $TEMP_DIR/test_pipeline_base/install/snpEff
[Parameters]
GENOME_NAME = test_genome
THREADS = 4
[Links]
FASTA_URL = http://example.com/mock.fasta.gz
GFF_URL = http://example.com/mock.gff.gz
SNPEFF_URL = http://example.com/mock_snpEff.zip
KRAKEN2_DB_URL = http://example.com/mock_k2.tar.gz
EOF

    # Create a mock Python script that reads the config file and returns values
    mkdir -p "$TEMP_DIR/scripts"
    cat <<'EOF' > "$TEMP_DIR/scripts/get_config_value.py"
#!/usr/bin/env python3
import sys
import configparser

def get_config_value(config_file, section, key):
    config = configparser.ConfigParser()
    config.read(config_file)
    return config.get(section, key) if config.has_option(section, key) else ''

if __name__ == '__main__':
    config_path = sys.argv[1]
    section_name = sys.argv[2]
    key_name = sys.argv[3]
    print(get_config_value(config_path, section_name, key_name))
EOF
    chmod +x "$TEMP_DIR/scripts/get_config_value.py"

    # Create mock executable commands that your script depends on
    mkdir -p "$TEMP_DIR/bin"

    # Mock wget and curl to simulate downloading files
    cat <<EOF > "$TEMP_DIR/bin/wget"
#!/bin/bash
touch "\$3"
echo "Mock wget called for \$1"
EOF
    cat <<EOF > "$TEMP_DIR/bin/curl"
#!/bin/bash
touch "\$3"
echo "Mock curl called for \$1"
EOF
    # Mock gunzip
    cat <<EOF > "$TEMP_DIR/bin/gunzip"
#!/bin/bash
# Check if the file exists before touching
if [ -f "\$2" ]; then
    filename="\$(basename "\$2")"
    unzipped_filename="\${filename%.gz}"
    touch "\$(dirname "\$2")/\$unzipped_filename"
    echo "Mock gunzip called for \$2"
    exit 0
fi
echo "Mock gunzip: File \$2 not found"
exit 1
EOF
    # Mock unzip
    cat <<EOF > "$TEMP_DIR/bin/unzip"
#!/bin/bash
# Simulate unzipping by creating a dummy jar file and a mock directory.
if [ -f "\$2" ]; then
    mkdir -p "\$4/snpEff"
    touch "\$4/snpEff/snpEff.jar"
    echo "Mock unzip called for \$2"
    exit 0
fi
echo "Mock unzip: File \$2 not found"
exit 1
EOF
    # Mock tar
    cat <<EOF > "$TEMP_DIR/bin/tar"
#!/bin/bash
if [ -f "\$2" ]; then
    mkdir -p "\$4/k2"
    echo "Mock tar called for \$2"
    exit 0
fi
echo "Mock tar: File \$2 not found"
exit 1
EOF
    # Mock bwa
    cat <<EOF > "$TEMP_DIR/bin/bwa"
#!/bin/bash
# Simulate bwa index creation by creating dummy index files.
if [ -f "\$2" ]; then
    touch "\$2.amb"
    touch "\$2.ann"
    touch "\$2.bwt"
    touch "\$2.sa"
    touch "\$2.pac"
    echo "Mock bwa index called for \$2"
    exit 0
fi
echo "Mock bwa: Reference file \$2 not found"
exit 1
EOF
    # Mock cnvpytor
    cat <<EOF > "$TEMP_DIR/bin/cnvpytor"
#!/bin/bash
echo "Mock cnvpytor -download called"
EOF
    # Mock java
    cat <<EOF > "$TEMP_DIR/bin/java"
#!/bin/bash
# Mock the snpEff build command.
echo "Mock java -jar snpEff build called"
# This mock no longer incorrectly tries to modify a config file.
EOF
    # Make all mock binaries executable
    chmod +x "$TEMP_DIR/bin/"*

    echo "Mock environment setup complete."
}

# ==============================================================================
# TEST RUNNER
# ==============================================================================

# Function to run the actual test
run_test() {
    echo "Running installation script with mock environment..."

    # Change to the temporary directory.
    cd "$TEMP_DIR" || { echo "Failed to change directory to $TEMP_DIR"; exit 1; }

    # Set up a fake PATH to prioritize our mock binaries.
    export PATH="$TEMP_DIR/bin:$PATH"

    # Run the installation script.
    # The 'source' command is used to run the script in the current shell,
    # allowing it to see the mock binaries on the PATH.
    if ! source "$OLDPWD/setup.sh"; then
        echo "FAIL: The installation script exited with an error."
        exit 1
    fi
    echo "Installation script finished."

    # Define the variables based on the mock config
    BASE_DIR="$TEMP_DIR/test_pipeline_base"
    DATA_DIR="$BASE_DIR/data"
    INSTALL_DIR="$BASE_DIR/install"
    BIN_DIR="$BASE_DIR/bin"
    SNPEFF_DIR="$INSTALL_DIR/snpEff"
    KRAKEN_DB="$INSTALL_DIR/kraken_db"
    GENOME_NAME="test_genome"

    # ==============================================================================
    # VERIFICATION CHECKS
    # ==============================================================================
    echo "Verifying installation..."
    
    # Check if all required directories exist
    echo "Checking for directories..."
    if [ ! -d "$BASE_DIR" ]; then echo "FAIL: BASE_DIR not created."; exit 1; fi
    if [ ! -d "$DATA_DIR" ]; then echo "FAIL: DATA_DIR not created."; exit 1; fi
    if [ ! -d "$INSTALL_DIR" ]; then echo "FAIL: INSTALL_DIR not created."; exit 1; fi
    if [ ! -d "$BIN_DIR" ]; then echo "FAIL: BIN_DIR not created."; exit 1; fi
    if [ ! -d "$SNPEFF_DIR" ]; then echo "FAIL: SNPEFF_DIR not created."; exit 1; fi
    if [ ! -d "$KRAKEN_DB" ]; then echo "FAIL: KRAKEN_DB not created."; exit 1; fi
    echo "All core directories exist."

    # Check for downloaded and processed reference genome files
    echo "Checking for reference genome files..."
    if [ ! -f "$DATA_DIR/reference.fasta" ]; then echo "FAIL: reference.fasta not created."; exit 1; fi
    if [ ! -f "$DATA_DIR/reference.gff" ]; then echo "FAIL: reference.gff not created."; exit 1; fi
    echo "Reference genome files exist."

    # Check for snpEff installation and custom genome files
    echo "Checking for snpEff files..."
    if [ ! -f "$SNPEFF_DIR/snpEff/snpEff.jar" ]; then echo "FAIL: snpEff.jar not created."; exit 1; fi
    if [ ! -d "$SNPEFF_DIR/snpEff/data/$GENOME_NAME" ]; then echo "FAIL: snpEff custom genome directory not created."; exit 1; fi
    if [ ! -f "$SNPEFF_DIR/snpEff/data/$GENOME_NAME/sequences.fa" ]; then echo "FAIL: snpEff sequences.fa not created."; exit 1; fi
    if [ ! -f "$SNPEFF_DIR/snpEff/data/$GENOME_NAME/genes.gff" ]; then echo "FAIL: snpEff genes.gff not created."; exit 1; fi
    # We can't easily test if the snpEff.config was modified, as we are mocking java
    echo "snpEff files exist."

    # Check for BWA index files
    echo "Checking for BWA index files..."
    if [ ! -f "$DATA_DIR/reference.fasta.amb" ]; then echo "FAIL: BWA index file '.amb' not created."; exit 1; fi
    if [ ! -f "$DATA_DIR/reference.fasta.ann" ]; then echo "FAIL: BWA index file '.ann' not created."; exit 1; fi
    # Add checks for other BWA index files if needed
    echo "BWA index files exist."

    # Check for Kraken2 database files
    echo "Checking for Kraken2 database files..."
    if [ ! -d "$KRAKEN_DB/k2" ]; then echo "FAIL: Kraken2 database directory not created."; exit 1; fi
    echo "Kraken2 database files exist."

    echo "All verification checks passed! Installation script works as expected."
    return 0
}

# ==============================================================================
# MAIN SCRIPT EXECUTION
# ==============================================================================

# Call the setup function
setup_mocks

# Call the test runner
run_test