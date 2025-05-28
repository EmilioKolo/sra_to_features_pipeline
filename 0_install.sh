#!/bin/bash

apt-get update
apt-get install -y bwa samtools bcftools sra-toolkit tabix bedtools

pip install pybedtools
pip install numpy==1.26.4
# Reinstall pandas
pip install --force-reinstall pandas==2.2.3

# Genome downloads and directory setup

# Create data directory
mkdir -p /content/data
# Create temporary directory
mkdir -p /content/data/tmp

# Download Reference Genome

# Define reference url folder and file
genome_url="https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gff_url="https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz"
# Define fasta reference and corresponding gff file
fasta_file="/content/data/reference.fasta"
gff_file="/content/data/reference.gff"

echo "Downloading reference genome..."
wget -O {fasta_file}.gz {genome_url}
wget -O {gff_file}.gz {gff_url}
gunzip -f {fasta_file}.gz
gunzip -f {gff_file}.gz
echo "Reference genome downloaded and unzipped to: {fasta_file}"

# Define snpEff variables
snpeff_dir="/content/data/bin"
genome_name="custom_ref"
snpeff_url='https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip'

# snpEff installation and creation of custom database
mkdir -p {snpeff_dir}

echo "Downloading and installing snpEff..."
wget -O {snpeff_dir}/snpEff.zip {snpeff_url}
unzip -o {snpeff_dir}/snpEff.zip -d {snpeff_dir}

echo "Creating snpEff custom genome..."
mkdir -p {snpeff_dir}/snpEff/data/{genome_name}
cp {fasta_file} {snpeff_dir}/snpEff/data/{genome_name}/sequences.fa
cp {gff_file} {snpeff_dir}/snpEff/data/{genome_name}/genes.gff
echo {genome_name}".genome : Custom genome" >> {snpeff_dir}/snpEff/snpEff.config
java -Xmx4g -jar {snpeff_dir}/snpEff/snpEff.jar build -gff3 -v {genome_name}
echo "snpEff custom genome generated."

# BWA indexing of the human genome
echo "Indexing reference genome: {fasta_file} with bwa..."
bwa index {fasta_file}
echo "Reference genome indexing complete."

# Generate genome.sizes file
python3 install.py

# Install cnvpytor
echo "Installing cnvpytor..."
git clone https://github.com/abyzovlab/CNVpytor.git /content/data
cd /content/data/CNVpytor; python setup.py install --user
echo "CNVpytor installed."
