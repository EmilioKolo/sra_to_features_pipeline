# SRA to features pipeline

Pipeline to go from an SRA ID (or fastq files) to features used to train an LLM.


# Pipeline schema

![Simplified pipeline ran in these scripts.](/pipeline_full.png)


# Instructions (Docker)

The code is optimised to run inside a Docker container, which handles all of the required programs and files. To change specific versions and names of folders, modify config.ini variables.
In order to run the pipeline, first create the container by running the following script:

`docker build -t features-pipeline .`

Then you can run the SRA IDs one by one with the following command:

`docker run -v {output_directory}:/content/data/output features-pipeline --sra-id {sra_id}`

All of the code is intended to work on cfDNA sequencing using Illumina with either single end or paired end sequencing reads.

Alternatively, the pipeline can be run with fastq files using the following command:

`docker run -v {output_directory}:/content/data/output -v {fastq_folder}:{fastq_folder} features-pipeline --fastq {fastq_file1} [{fastq_file2}]`

The variables `fastq_file1` and `fastq_file2` must be the full name of fastq files (including the directory), and they both must be in the directory `fastq_folder`.

Using only `fastq_file1` assumes single-end Illumina sequencing. Defining `fastq_file2` assumes paired-end Illumina sequencing.


# Instructions (no Docker)

To run the code without a Docker container, first it is recommended to define the BASE_DIR value in config.ini to a local folder, then run the following command from inside this folder to install the required programs:

`sudo chmod +x setup.sh; sudo setup.sh; python3 install.py`

After the setup finishes preparing the needed programs and files, the pipeline can be run with the following command:

`python3 main.py --sra-id {sra_id} --output-dir {output_dir} --threads {threads}`

Where `sra_id` corresponds to a SRA ID, `output_dir` to a local folder where the output files will be created, and `threads` to the number of threads to be used by some functions with multithreading.

Alternatively, the pipeline can be run with fastq files using the following command:

`python3 main.py --fastq {fastq_file1} [{fastq_file2}] --output-dir {output_dir} --threads {threads}`

The variables `fastq_file1` and `fastq_file2` must be the full name of fastq files, and both of them need to be in the same folder.

Using only `fastq_file1` assumes single-end Illumina sequencing. Defining `fastq_file2` assumes paired-end Illumina sequencing.


# Output

The script will generate a file with the following features:
* Fragment lengths and their mean, median and st. deviation (for paired-end sequencing)
* Number of variants in genome bins through the entire genome
* Number of variants per region/gene in selected regions and all genes
* Synonymous/Nonsynonymous variant proportion per gene for all genes
* Number of CNVs per chromosome

It also outputs vcf files and quality control results.

