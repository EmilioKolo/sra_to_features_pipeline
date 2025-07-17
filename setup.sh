#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# Define the config file
CONFIG_FILE="config.ini"
# Get relevant variables
BASE_DIR=$(python3 get_config_value.py "$CONFIG_FILE" "Paths" "BASE_DIR")
if [ -z "$BASE_DIR" ]; then
    echo "Error: BASE_DIR not found in config.ini"
    exit 1
fi
THREADS=$(python3 get_config_value.py "$CONFIG_FILE" "Parameters" "THREADS")
if [ -z "$THREADS" ]; then
    echo "Error: THREADS not found in config.ini"
    exit 1
fi
GENOME_NAME=$(python3 get_config_value.py "$CONFIG_FILE" "Parameters" "GENOME_NAME")
if [ -z "$GENOME_NAME" ]; then
    echo "Error: GENOME_NAME not found in config.ini"
    exit 1
fi
FASTA_URL=$(python3 get_config_value.py "$CONFIG_FILE" "Links" "FASTA_URL")
if [ -z "$FASTA_URL" ]; then
    echo "Error: FASTA_URL not found in config.ini"
    exit 1
fi
GFF_URL=$(python3 get_config_value.py "$CONFIG_FILE" "Links" "GFF_URL")
if [ -z "$GFF_URL" ]; then
    echo "Error: GFF_URL not found in config.ini"
    exit 1
fi
SNPEFF_URL=$(python3 get_config_value.py "$CONFIG_FILE" "Links" "SNPEFF_URL")
if [ -z "$SNPEFF_URL" ]; then
    echo "Error: SNPEFF_URL not found in config.ini"
    exit 1
fi
KRAKEN2_DB_URL=$(python3 get_config_value.py "$CONFIG_FILE" "Links" "KRAKEN2_DB_URL")
if [ -z "$KRAKEN2_DB_URL" ]; then
    echo "Error: KRAKEN2_DB_URL not found in config.ini"
    exit 1
fi

# Make sure that BASE_DIR exists
mkdir -p "$BASE_DIR"

# Create data, install, logs, bin and tmp directories
DATA_DIR="$BASE_DIR/data"
INSTALL_DIR="$BASE_DIR/install"
LOGS_DIR="$DATA_DIR/logs"
BIN_DIR="$DATA_DIR/bin"
TMP_DIR="$DATA_DIR/tmp"
mkdir -p "$DATA_DIR"
mkdir -p "$INSTALL_DIR"
mkdir -p "$LOGS_DIR"
mkdir -p "$BIN_DIR"
mkdir -p "$TMP_DIR"
# Add the bin directory to PATH
export PATH="$BIN_DIR:$PATH"

# Function to download with either curl or wget (whichever exists)
download() {
    if command -v curl >/dev/null; then
        curl -L "$1" -o "$2"
    elif command -v wget >/dev/null; then
        wget "$1" -O "$2"
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
# Define snpEff variables
snpeff_dir="$BIN_DIR"
genome_name="$GENOME_NAME"
snpeff_jar="$snpeff_dir/snpEff/snpEff.jar"
# Install snpEff
mkdir -p "$snpeff_dir"
download "$SNPEFF_URL" "$snpeff_dir/snpEff.zip"
unzip -o "$snpeff_dir/snpEff.zip" -d "$snpeff_dir"
rm -f "$snpeff_dir/snpEff.zip"
log "Finished snpEff download. Starting snpEff custom genome setup."
# Create custom genome for snpeff
mkdir -p "$snpeff_dir/snpEff/data/$genome_name"
cp "$DATA_DIR/reference.fasta" "$snpeff_dir/snpEff/data/$genome_name/sequences.fa"
cp "$DATA_DIR/reference.gff" "$snpeff_dir/snpEff/data/$genome_name/genes.gff"
echo "${genome_name}.genome : Custom genome" >> "$snpeff_dir/snpEff/snpEff.config"
java -Xmx4g -jar "$snpeff_jar" build -gff3 -v "$genome_name"

log "Finished snpEff custom genome setup. Starting BWA indexing."
# Index reference genome with bwa
bwa index "$DATA_DIR/reference.fasta"

log "Finished BWA indexing. Starting cnvpytor data download."
# Download cnvpytor data
cnvpytor -download

log "Finished cnvpytor data download. Starting Kraken2 database creation."
# Download Kraken2 database
kraken_db="$INSTALL_DIR/kraken2-db"
mkdir -p "$kraken_db"
download "$KRAKEN2_DB_URL" "$kraken_db/k2.tar.gz"
tar -xzf "$kraken_db/k2.tar.gz" -C "$kraken_db"
rm -f "$kraken_db/k2.tar.gz"
log "Finished Kraken2 database creation."
