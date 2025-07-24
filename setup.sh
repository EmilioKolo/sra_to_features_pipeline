#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# Define the config file
CONFIG_FILE="config.ini"
# Get relevant variables
BASE_DIR=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Paths" "BASE_DIR")
if [ -z "$BASE_DIR" ]; then
    echo "Error: BASE_DIR not found in config.ini"
    exit 1
fi
KRAKEN_DB=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Paths" "KRAKEN_DB")
if [ -z "$KRAKEN_DB" ]; then
    echo "Error: KRAKEN_DB not found in config.ini"
    exit 1
fi
SNPEFF_DIR=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Paths" "SNPEFF_DIR")
if [ -z "$SNPEFF_DIR" ]; then
    echo "Error: SNPEFF_DIR not found in config.ini"
    exit 1
fi
GENOME_NAME=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Parameters" "GENOME_NAME")
if [ -z "$GENOME_NAME" ]; then
    echo "Error: GENOME_NAME not found in config.ini"
    exit 1
fi
THREADS=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Parameters" "THREADS")
if [ -z "$THREADS" ]; then
    echo "Error: THREADS not found in config.ini"
    exit 1
fi
FASTA_URL=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Links" "FASTA_URL")
if [ -z "$FASTA_URL" ]; then
    echo "Error: FASTA_URL not found in config.ini"
    exit 1
fi
GFF_URL=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Links" "GFF_URL")
if [ -z "$GFF_URL" ]; then
    echo "Error: GFF_URL not found in config.ini"
    exit 1
fi
SNPEFF_URL=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Links" "SNPEFF_URL")
if [ -z "$SNPEFF_URL" ]; then
    echo "Error: SNPEFF_URL not found in config.ini"
    exit 1
fi
KRAKEN2_DB_URL=$(python3 scripts/get_config_value.py "$CONFIG_FILE" "Links" "KRAKEN2_DB_URL")
if [ -z "$KRAKEN2_DB_URL" ]; then
    echo "Error: KRAKEN2_DB_URL not found in config.ini"
    exit 1
fi

# Make sure that BASE_DIR exists
mkdir -p "$BASE_DIR"

# Create data, install, logs, bin and tmp directories
DATA_DIR="$BASE_DIR/data"
INSTALL_DIR="$BASE_DIR/install"
BIN_DIR="$BASE_DIR/bin"
mkdir -p "$DATA_DIR"
mkdir -p "$INSTALL_DIR"
mkdir -p "$BIN_DIR"
# Define snpEff variables
snpeff_jar="$SNPEFF_DIR/snpEff/snpEff.jar"
# Add the bin directory to PATH
export PATH="$BIN_DIR:$PATH"

# Function to download with either curl or wget (whichever exists)
download() {
    if command -v wget >/dev/null; then
        wget "$1" -O "$2"
    elif command -v curl >/dev/null; then
        curl -L "$1" -o "$2"
    else
        echo "ERROR: Need either curl or wget to download files."
        exit 1
    fi
}
# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log "Starting reference genome download."
# Download reference genome
download "$FASTA_URL" "$DATA_DIR/reference.fasta.gz"
gunzip -f "$DATA_DIR/reference.fasta.gz"
download "$GFF_URL" "$DATA_DIR/reference.gff.gz"
gunzip -f "$DATA_DIR/reference.gff.gz"

log "Finished reference genome download. Starting snpEff download."
# Install snpEff
mkdir -p "$SNPEFF_DIR"
download "$SNPEFF_URL" "$SNPEFF_DIR/snpEff.zip"
unzip -o "$SNPEFF_DIR/snpEff.zip" -d "$SNPEFF_DIR"
rm -f "$SNPEFF_DIR/snpEff.zip"
log "Finished snpEff download. Starting snpEff custom genome setup."
# Create custom genome for snpeff
mkdir -p "$SNPEFF_DIR/snpEff/data/$GENOME_NAME"
cp "$DATA_DIR/reference.fasta" "$SNPEFF_DIR/snpEff/data/$GENOME_NAME/sequences.fa"
cp "$DATA_DIR/reference.gff" "$SNPEFF_DIR/snpEff/data/$GENOME_NAME/genes.gff"
echo "${GENOME_NAME}.genome : Custom genome" >> "$SNPEFF_DIR/snpEff/snpEff.config"
java -Xmx4g -jar "$snpeff_jar" build -gff3 -v "$GENOME_NAME"

log "Finished snpEff custom genome setup. Starting BWA indexing."
# Index reference genome with bwa
bwa index "$DATA_DIR/reference.fasta"

log "Finished BWA indexing. Starting cnvpytor data download."
# Download cnvpytor data
cnvpytor -download

log "Finished cnvpytor data download. Starting Kraken2 database creation."
# Download Kraken2 database
mkdir -p "$KRAKEN_DB"
download "$KRAKEN2_DB_URL" "$KRAKEN_DB/k2.tar.gz"
tar -xzf "$KRAKEN_DB/k2.tar.gz" -C "$KRAKEN_DB"
rm -f "$KRAKEN_DB/k2.tar.gz"
log "Finished Kraken2 database creation."
