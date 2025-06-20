#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# Get environment variables
source config.env

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

# Download reference genome
download "$FASTA_URL" "$DATA_DIR/reference.fasta.gz"
gunzip -f "$DATA_DIR/reference.fasta.gz"
download "$GFF_URL" "$DATA_DIR/reference.gff.gz"
gunzip -f "$DATA_DIR/reference.gff.gz"

# Define snpEff variables
snpeff_dir="$BIN_DIR"
genome_name="$GENOME_NAME"
snpeff_jar="$snpeff_dir/snpEff/snpEff.jar"
# Install snpEff
mkdir -p "$snpeff_dir"
download "$SNPEFF_URL" "$snpeff_dir/snpEff.zip"
unzip -o "$snpeff_dir/snpEff.zip" -d "$snpeff_dir"
rm -f "$snpeff_dir/snpEff.zip"
# Create custom genome for snpeff
mkdir -p "$snpeff_dir/snpEff/data/$genome_name"
cp "$DATA_DIR/reference.fasta" "$snpeff_dir/snpEff/data/$genome_name/sequences.fa"
cp "$DATA_DIR/reference.gff" "$snpeff_dir/snpEff/data/$genome_name/genes.gff"
echo "${genome_name}.genome : Custom genome" >> "$snpeff_dir/snpEff/snpEff.config"
java -Xmx4g -jar "$snpeff_jar" build -gff3 -v "$genome_name"

# Index reference genome with bwa
bwa index "$DATA_DIR/reference.fasta"

# Download cnvpytor data
cnvpytor -download

# Download Kraken2 database
kraken_db="$INSTALL_DIR/kraken2-db"
mkdir -p "$kraken_db"
download "$KRAKEN2_DB_URL" "$kraken_db/k2.tar.gz"
tar -xzf "$kraken_db/k2.tar.gz" -C "$kraken_db"
rm -f "$kraken_db/k2.tar.gz"
