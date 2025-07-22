#!/usr/bin/env python3


import os


### Modular function ###

def perform_checks(dict_var:dict) -> None:
    """
    Performs checks on the different intermediate files.
    """
    # Define a general output folder for checks
    checks_out = os.path.join(dict_var['OUTPUT_DIR'] + 'checks')
    # Make the output directory if it does not exist
    try:
        os.mkdir(checks_out)
    except FileExistsError:
        pass

    # Define other needed output folders
    kraken_out = os.path.join(checks_out, 'kraken_out')
    fastqc_out = os.path.join(checks_out, 'fastqc_out')

    try:
        print("Classifying unaligned reads with Kraken2...")
        # Perform Kraken check on sorted BAM file
        check_kraken(
            dict_var['sra_id'],
            bam_file=dict_var['sorted_bam'],
            kraken_db=dict_var['kraken_db'],
            output_dir=kraken_out,
            threads=dict_var['THREADS']
        )
    except Exception as e:
        print(f"Error classifying unaligned reads: {e}")
        print("Skipping unaligned reads classification.")
    
    try:
        print('Analyzing fastq files with FastQC')
        # Go through full paths of all fastq files
        for fastq_file in dict_var['l_fastq_full']:
            # Run FastQC on each file
            run_fastqc(fastq_file, fastqc_out)
    except Exception as e:
        print(f"Error analyzing fastq files: {e}")
        print("Skipping fastq file analysis.")
    return None


### Shorter functions ###

def check_kraken(
        sra_id:str,
        bam_file:str,
        kraken_db:str,
        output_dir:str,
        threads:int=1
    ) -> None:
    """
    Classify unaligned reads from a BAM file using Kraken2.
    """
    # Make the output directory if it does not exist
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    # Define intermediate and output files
    unaligned_bam = os.path.join(output_dir, sra_id+'_unaligned.bam')
    unaligned_fastq = os.path.join(output_dir, sra_id+'_unaligned.fastq')
    kraken_output = os.path.join(output_dir, sra_id+'_kraken2_output.txt')
    kraken_report = os.path.join(output_dir, sra_id+'_kraken2_report.txt')

    # Extract unaligned reads
    print('Extracting unaligned reads from BAM file...')
    l = f'samtools view -f 4 -b {bam_file} -o {unaligned_bam}'
    os.system(l)
    print('Converting unaligned BAM to FASTQ...')
    l = f'samtools fastq {unaligned_bam} -o {unaligned_fastq}'
    os.system(l)

    print('Running Kraken2 classification...')
    # Classify with Kraken2
    l = f'kraken2 --db {kraken_db} --threads {threads}'
    l += f' --output {kraken_output} --report {kraken_report}'
    l += f' {unaligned_fastq}'
    os.system(l)

    # Visualize output data
    p = 'Kraken2 classification complete. Results in:'
    p += f'\n * {kraken_output}\n * {kraken_report}'
    print(p)

    # Clean up intermediate files
    os.remove(unaligned_bam)
    os.remove(unaligned_fastq)
    return None

def run_fastqc(
        fastq_file:str,
        output_dir:str
    ) -> None:
    """
    Runs FastQC on a fastq file.
    """
    # Make the output directory if it does not exist
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    # Initialize script
    l = 'fastqc -q'
    # Define output dir
    l += f' -o {output_dir}'
    # Define fastq file to analyze
    l += f' {fastq_file}'
    # Run the script
    os.system(l)
    return None

def run_fastqc_folder(
        fastq_folder:str,
        output_dir:str,
        threads:int=1
    ) -> None:
    """
    Runs FastQC on a folder with fastq files.
    """
    # Make the output directory if it does not exist
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    # Initialize script
    l = 'fastqc -q'
    # Define output dir
    l += f' -o {output_dir}'
    # Define thread number
    l += f' -t {threads}'
    # Define fastq folder to analyze
    l += f' {fastq_folder}/*'
    # Run the script
    os.system(l)
    return None
