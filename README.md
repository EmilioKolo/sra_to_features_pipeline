# SRA to features pipeline
Pipeline to go from an SRA ID to features used to train an LLM.

# Instructions

- Run 0_install.sh as executable (with super user permissions). It downloads the Homo_sapiens.GRCh38.109 genome from Ensembl by default.
- Run 1_main.py with python passing a valid SRA ID with either single end or paired end sequencing reads. Intended to work on cfDNA sequencing using Illumina.

# Output

The script will generate a file with the following features:
* Small fragment lengths and their mean, median and st. deviation
* Number of variants in genome bins through the entire genome
* Number of variants per region/gene in selected regions and all genes
* Synonymous/Nonsynonymous variant proportion per gene for all genes
* Number of CNVs per chromosome
