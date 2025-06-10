# SRA to features pipeline
Pipeline to go from an SRA ID to features used to train an LLM.

# Instructions

One-line run command (runs with conda):
- Run "chmod +x ./0_main.sh; ./0_main.sh". It makes 0_main.sh an executable and runs it, which sets up the conda environment and runs 2_run_multiple.py.

Granular run commands (needs super user permissions and runs without conda):
- Run 0_install.sh as executable (with super user permissions). It downloads the Homo_sapiens.GRCh38.109 genome from Ensembl by default.
- Run 1_main.py or 2_run_multiple.py. For 1_main.py, you need to pass a valid SRA ID as a system variable. 2_run_multiple.py takes a file path to a list of SRA IDs; if none is provided, it works with sra_table_selected.txt in the Git folder.

All of the code is intended to work on cfDNA sequencing using Illumina with either single end or paired end sequencing reads.

# Output

The script will generate a file with the following features:
* Small fragment lengths and their mean, median and st. deviation
* Number of variants in genome bins through the entire genome
* Number of variants per region/gene in selected regions and all genes
* Synonymous/Nonsynonymous variant proportion per gene for all genes
* Number of CNVs per chromosome
