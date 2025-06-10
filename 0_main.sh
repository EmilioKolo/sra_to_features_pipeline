#!/bin/bash

# Exit on error and undefined variables
set -euo pipefail

# Define and make install_logs folder
install_logs="$HOME/install_logs"
mkdir -p "$install_logs"

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Check and install Miniconda if not already installed
if [[ ! -d "$HOME/miniconda" ]]; then
    log "Downloading and installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-py310_25.3.1-1-Linux-x86_64.sh -O "$HOME/miniconda.sh"
    bash "$HOME/miniconda.sh" -b -p "$HOME/miniconda" > "$install_logs/miniconda_bash.log" 2>&1
    rm -f "$HOME/miniconda.sh"
else
    log "Miniconda already installed at $HOME/miniconda, skipping installation."
fi
# Initialize conda
source "$HOME/miniconda/bin/activate"

# Check and create conda environment if not exists
if ! conda env list | grep -q 'sra_to_feats'; then
    log "Creating conda environment sra_to_feats..."
    conda create -n sra_to_feats python=3.10 -y > "$install_logs/conda_create.log" 2>&1
else
    log "Conda environment sra_to_feats already exists, skipping creation."
fi

# Install conda packages if environment was just created or packages might be missing
if [[ ! -d "$HOME/miniconda/envs/sra_to_feats" ]] || \
   [[ ! $(conda list -n sra_to_feats | grep bcftools) ]]; then
    log "Installing conda packages..."
    conda install -n sra_to_feats -c conda-forge \
        -y bioconda::bcftools \
        bioconda::bedtools \
        bioconda::bwa \
        bioconda::kraken2 \
        bioconda::samtools \
        bioconda::sra-tools=3.0.0 \
        bioconda::tabix \
        conda-forge::unzip \
        cyclus::java-jre > "$install_logs/conda_install.log" 2>&1
else
    log "Conda packages already installed, skipping."
fi

# Activate environment
conda activate sra_to_feats > "$install_logs/conda_activate.log" 2>&1

# Install Python packages if not already installed
log "Checking Python packages..."
python3 -m pip install --upgrade pip
python3 -m pip install numpy==1.26.4 pandas==2.1.4 pybedtools requests cnvpytor==1.3.1

# Create data and logging directories
mkdir -p "$HOME/content/data"
mkdir -p "$HOME/content/data/logs"

# Download reference genome if it does not exist
if [[ ! -f "$HOME/content/data/reference.fasta" ]]; then
    log "Downloading reference genome..."
    wget -O "$HOME/content/data/reference.fasta.gz" https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip -f "$HOME/content/data/reference.fasta.gz"
else
    log "Reference genome already exists, skipping download."
fi
if [[ ! -f "$HOME/content/data/reference.gff" ]]; then
    log "Downloading reference gff..."
    wget -O "$HOME/content/data/reference.gff.gz" https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz
    gunzip -f "$HOME/content/data/reference.gff.gz"
else
    log "Reference gff already exists, skipping download."
fi

# Install snpEff if it does not exist
snpeff_dir="$HOME/content/data/bin"
genome_name="custom_ref"
snpeff_jar="$snpeff_dir/snpEff/snpEff.jar"

if [[ ! -f "$snpeff_jar" ]]; then
    log "Downloading and installing snpEff..."
    mkdir -p "$snpeff_dir"
    wget -O "$snpeff_dir/snpEff.zip" https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
    unzip -o "$snpeff_dir/snpEff.zip" -d "$snpeff_dir"
    rm -f "$snpeff_dir/snpEff.zip"
else
    log "snpEff already installed, skipping download."
fi

# Create custom genome if it does not exist
if [[ ! -d "$snpeff_dir/snpEff/data/$genome_name" ]]; then
    log "Creating snpEff custom genome..."
    mkdir -p "$snpeff_dir/snpEff/data/$genome_name"
    cp "$HOME/content/data/reference.fasta" "$snpeff_dir/snpEff/data/$genome_name/sequences.fa"
    cp "$HOME/content/data/reference.gff" "$snpeff_dir/snpEff/data/$genome_name/genes.gff"
    echo "${genome_name}.genome : Custom genome" >> "$snpeff_dir/snpEff/snpEff.config"
    java -Xmx4g -jar "$snpeff_jar" build -gff3 -v "$genome_name" > "$install_logs/snpeff.log" 2>&1
else
    log "snpEff custom genome already exists, skipping creation."
fi

# Index reference genome if not indexed
if [[ ! -f "$HOME/content/data/reference.fasta.bwt" ]]; then
    log "Indexing reference genome with bwa..."
    bwa index "$HOME/content/data/reference.fasta" > "$install_logs/bwa_index.log" 2>&1
else
    log "Reference genome already indexed, skipping."
fi

# Run install.py if output does not exist
if [[ ! -f "$HOME/content/data/genome.sizes" ]]; then
    log "Running install.py..."
    python3 install.py "$HOME" > "$install_logs/install_py.log" 2>&1
else
    log "install.py was already run, skipping."
fi

# Download cnvpytor data if it does not exist
if [[ ! -d "$HOME/.cnvpytor/data" ]]; then
    log "Downloading cnvpytor data..."
    cnvpytor -download > "$install_logs/cnvpytor.log" 2>&1
else
    log "CNVpytor data already exists, skipping download."
fi

# Download Kraken2 database if it does not exist
kraken_db="$HOME/kraken2-db"
if [[ ! -f "$kraken_db/hash.k2d" ]]; then
    log "Downloading Kraken2 database..."
    mkdir -p "$kraken_db"
    wget -O "$kraken_db/k2.tar.gz" https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz
    tar -xvzf "$kraken_db/k2.tar.gz" -C "$kraken_db"
    rm -f "$kraken_db/k2.tar.gz"
else
    log "Kraken2 database already exists, skipping download."
fi

# Run the main script
log "Starting main script execution..."
python3 2_run_multiple.py
