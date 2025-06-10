#!/bin/bash


# Define and make install_logs folder
install_logs=$HOME/install_logs
mkdir -p $install_logs

# Download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_25.3.1-1-Linux-x86_64.sh -O "$HOME/miniconda.sh"
bash "$HOME/miniconda.sh" -b -p "$HOME/miniconda" > $install_logs/miniconda_bash.log 2>&1
source "$HOME/miniconda/bin/activate"

# Create environment and install packages
conda create -n sra_to_feats python=3.10 -y > $install_logs/conda_create.log 2>&1
conda install -n sra_to_feats -c conda-forge \
  -y bioconda::bcftools \
  bioconda::bedtools \
  bioconda::bwa \
  bioconda::kraken2 \
  bioconda::samtools \
  bioconda::sra-tools=3.0.0 \
  bioconda::tabix \
  conda-forge::unzip \
  cyclus::java-jre > $install_logs/conda_install.log 2>&1
conda activate sra_to_feats > $install_logs/conda_activate.log 2>&1

# Python packages
python3 -m pip install numpy==1.26.4
python3 -m pip install pandas==2.1.4
python3 -m pip install pybedtools
python3 -m pip install requests
python3 -m pip install cnvpytor==1.3.1

# Genome downloads and directory setup

# Create data directory
mkdir -p "$HOME/content/data"
# Create logging directory
mkdir -p "$HOME/content/data/logs"

# Download Reference Genome

# Define reference url folder and file
genome_url="https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gff_url="https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz"
# Define fasta reference and corresponding gff file
fasta_file="$HOME/content/data/reference.fasta"
gff_file="$HOME/content/data/reference.gff"

echo "Downloading reference genome..."
wget -O "$HOME/content/data/reference.fasta.gz" https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -O "$HOME/content/data/reference.gff.gz" https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz
gunzip -f "$HOME/content/data/reference.fasta.gz"
gunzip -f "$HOME/content/data/reference.gff.gz"
echo "Reference genome downloaded and unzipped to: $HOME/content/data/reference.fasta"

# Define snpEff variables
snpeff_dir="$HOME/content/data/bin"
genome_name="custom_ref"
snpeff_url='https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip'

# snpEff installation and creation of custom database
mkdir -p "$HOME/content/data/bin"

echo "Downloading and installing snpEff..."
wget -O "$HOME/content/data/bin/snpEff.zip" https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
unzip -o "$HOME/content/data/bin/snpEff.zip" -d "$HOME/content/data/bin"

echo "Creating snpEff custom genome..."
mkdir -p "$HOME/content/data/bin/snpEff/data/$genome_name"
cp "$HOME/content/data/reference.fasta" "$HOME/content/data/bin/snpEff/data/$genome_name/sequences.fa"
cp "$HOME/content/data/reference.gff" "$HOME/content/data/bin/snpEff/data/$genome_name/genes.gff"
echo "${genome_name}.genome : Custom genome" >> "$HOME/content/data/bin/snpEff/snpEff.config"
java -Xmx4g -jar "$HOME/content/data/bin/snpEff/snpEff.jar build" -gff3 -v $genome_name > $install_logs/snpeff.log 2>&1
echo "snpEff custom genome generated."

# BWA indexing of the human genome
echo "Indexing reference genome: $HOME/content/data/reference.fasta with bwa..."
bwa index "$HOME/content/data/reference.fasta" > $install_logs/bwa_index.log 2>&1
echo "Reference genome indexing complete."

echo "Running install.py"
# Generate genome.sizes file
python3 install.py "$HOME" > $install_logs/install_py.log 2>&1

# Install cnvpytor
echo "Downloading cnvpytor data..."
cnvpytor -download > $install_logs/cnvpytor.log 2>&1
echo "CNVpytor installed."

# Download standard Kraken2 database (can take $HOME30-50 GB)
rm -rf "$HOME/kraken2-db"
mkdir -p "$HOME/kraken2-db"
wget -O "$HOME/kraken2-db/k2.tar.gz" https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz
tar -xvzf "$HOME/kraken2-db/k2.tar.gz"
rm "$HOME/kraken2-db/k2.tar.gz"
echo "Kraken2 database downloaded and extracted to $HOME/kraken2-db."

# Run the main script
python3 2_run_multiple.py
