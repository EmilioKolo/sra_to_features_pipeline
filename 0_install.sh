#!/bin/bash

# Update packages
sudo apt-get update
# Install system dependencies
sudo apt-get install -y \
  bwa \
  samtools \
  bcftools \
  sra-toolkit \
  tabix \
  bedtools \
  python3-pip \
  unzip \
  default-jre

# Python packages
python3 -m pip install numpy==1.26.4
python3 -m pip install pandas==2.1.4
python3 -m pip install pybedtools

# Genome downloads and directory setup

# Create data directory
mkdir -p ~/content/data
# Create temporary directory
mkdir -p ~/content/data/tmp

# Download Reference Genome

# Define reference url folder and file
genome_url="https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gff_url="https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz"
# Define fasta reference and corresponding gff file
fasta_file="~/content/data/reference.fasta"
gff_file="~/content/data/reference.gff"

echo "Downloading reference genome..."
wget -O ~/content/data/reference.fasta.gz https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -O ~/content/data/reference.gff.gz https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz
gunzip -f ~/content/data/reference.fasta.gz
gunzip -f ~/content/data/reference.gff.gz
echo "Reference genome downloaded and unzipped to: ~/content/data/reference.fasta"

# Define snpEff variables
snpeff_dir="~/content/data/bin"
genome_name="custom_ref"
snpeff_url='https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip'

# snpEff installation and creation of custom database
mkdir -p ~/content/data/bin

echo "Downloading and installing snpEff..."
wget -O ~/content/data/bin/snpEff.zip https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
unzip -o ~/content/data/bin/snpEff.zip -d ~/content/data/bin

echo "Creating snpEff custom genome..."
mkdir -p ~/content/data/bin/snpEff/data/genome_name
cp ~/content/data/reference.fasta ~/content/data/bin/snpEff/data/genome_name/sequences.fa
cp ~/content/data/reference.gff ~/content/data/bin/snpEff/data/genome_name/genes.gff
echo genome_name".genome : Custom genome" >> ~/content/data/bin/snpEff/snpEff.config
java -Xmx4g -jar ~/content/data/bin/snpEff/snpEff.jar build -gff3 -v genome_name
echo "snpEff custom genome generated."

# BWA indexing of the human genome
echo "Indexing reference genome: ~/content/data/reference.fasta with bwa..."
bwa index ~/content/data/reference.fasta
echo "Reference genome indexing complete."

# Generate genome.sizes file
python3 install.py

# Install cnvpytor
echo "Installing cnvpytor..."
git clone https://github.com/abyzovlab/CNVpytor.git ~/content/data
cd ~/content/data/CNVpytor; python setup.py install --user
echo "CNVpytor installed."
