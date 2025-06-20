# SRA to features pipeline
Pipeline to go from an SRA ID to features used to train an LLM.

# Instructions

The code is intended to run inside a Docker container, which handles all of the required programs and files. To change specific versions and names of folders, modify config.env variables.
In order to run the pipeline, first create the container by running the following script:

`$ docker build -t features-pipeline .`

Then you can run main.py and pass the path to a list of SRA IDs to analyze. Example script:

`$ python3 main.py ./sra_table_selected.txt {output_directory}`

The outputs will be generated in the provided output directory. If an output directory is not provided, the results will be on a folder next to the extracted git code.

OPTIONAL: Run the SRA IDs one by one with the following command:

`$ docker run -v {output_directory}:/content/data features-pipeline {sra_id} /content/data`

All of the code is intended to work on cfDNA sequencing using Illumina with either single end or paired end sequencing reads.

# Output

The script will generate a file with the following features:
* Small fragment lengths and their mean, median and st. deviation
* Number of variants in genome bins through the entire genome
* Number of variants per region/gene in selected regions and all genes
* Synonymous/Nonsynonymous variant proportion per gene for all genes
* Number of CNVs per chromosome
