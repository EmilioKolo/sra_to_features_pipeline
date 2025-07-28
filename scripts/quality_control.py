#!/usr/bin/env python3


"""
Functions used for quality control and quality checks of the output files.
"""

from scripts.log_scripts import *
import os


### Modular function ###

def perform_checks(dict_var:dict) -> None:
    """
    Performs checks on the different intermediate files.
    """
    # Define a general output folder for checks
    checks_out = os.path.join(dict_var['OUTPUT_DIR'], 'checks')
    # Make the output directory if it does not exist
    try:
        os.mkdir(checks_out)
    except FileExistsError:
        pass

    # Define other needed output folders
    kraken_out = os.path.join(checks_out, 'kraken_out')
    fastqc_out = os.path.join(checks_out, 'fastqc_out')
    coverage_out = os.path.join(checks_out, 'coverage_out')

    # FastQC analysis
    try:
        log_print(
            'Analyzing fastq files with FastQC...',
            level='info',
            log_file=dict_var['log_print']
        )
        # Go through full paths of all fastq files
        for fastq_file in dict_var['l_fastq_full']:
            # Run FastQC on each file
            run_fastqc(
                fastq_file,
                fastqc_out,
                log_scr=dict_var['log_scripts']
            )
    except Exception as e:
        log_print(
            f"Error analyzing fastq files: {e}",
            level='error',
            log_file=dict_var['log_print']
        )
        log_print(
            'Skipping fastq file analysis.',
            level='info',
            log_file=dict_var['log_print']
        )

    # Coverage analysis
    try:
        log_print(
            'Obtaining coverage statistics...',
            level='info',
            log_file=dict_var['log_print']
        )
        # Get coverage statistics using bedtools
        run_bedtools_coverage(
            dict_var['sra_id'],
            dict_var['sorted_bam'],
            coverage_out,
            log_scr=dict_var['log_scripts'],
            region_file=dict_var['bed_genes']
        )
        # Get coverage statistics using samtools
        run_samtools_depth(
            dict_var['sra_id'],
            dict_var['sorted_bam'],
            coverage_out,
            log_scr=dict_var['log_scripts'],
            region_file=dict_var['bed_genes']
        )
    except Exception as e:
        log_print(
            f"Error obtaining coverage statistics: {e}",
            level='error',
            log_file=dict_var['log_print']
        )
        log_print(
            'Skipping coverage statistics.',
            level='info',
            log_file=dict_var['log_print']
        )

    # Kraken2 analysis
    try:
        log_print(
            'Classifying unaligned reads with Kraken2...',
            level='info',
            log_file=dict_var['log_print']
        )
        # Perform Kraken check on sorted BAM file
        check_kraken(
            dict_var['sra_id'],
            bam_file=dict_var['sorted_bam'],
            kraken_db=dict_var['kraken_db'],
            output_dir=kraken_out,
            log_file=dict_var['log_print'],
            log_scr=dict_var['log_scripts'],
            threads=dict_var['THREADS']
        )
    except Exception as e:
        log_print(
            f"Error classifying unaligned reads: {e}",
            level='error',
            log_file=dict_var['log_print']
        )
        log_print(
            'Skipping unaligned reads classification.',
            level='info',
            log_file=dict_var['log_print']
        )
    
    return None


### Shorter functions ###

def check_kraken(
        sra_id:str,
        bam_file:str,
        kraken_db:str,
        output_dir:str,
        log_file:str,
        log_scr:str,
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
    log_print(
        'Extracting unaligned reads from BAM file...',
        level='info',
        log_file=log_file
    )
    l = f'samtools view -f 4 -b {bam_file} -o {unaligned_bam}'
    log_code(l, log_file=log_scr)
    os.system(l)
    log_print(
        'Converting unaligned BAM to FASTQ...',
        level='info',
        log_file=log_file
    )
    l = f'samtools fastq {unaligned_bam} -o {unaligned_fastq}'
    log_code(l, log_file=log_scr)
    os.system(l)

    log_print(
        'Running Kraken2 classification...',
        level='info',
        log_file=log_file
    )
    # Classify with Kraken2
    l = f'kraken2 --db {kraken_db} --threads {threads}'
    l += f' --output {kraken_output} --report {kraken_report}'
    l += f' {unaligned_fastq}'
    log_code(l, log_file=log_scr)
    os.system(l)

    # Visualize output data
    p = 'Kraken2 classification complete. Results in:'
    p += f'\n * {kraken_output}\n * {kraken_report}'
    log_print(
        p,
        level='info',
        log_file=log_file
    )

    # Clean up intermediate files
    os.remove(unaligned_bam)
    os.remove(unaligned_fastq)
    return None

def run_bedtools_coverage(
        sra_id:str,
        bam_sorted:str,
        output_dir:str,
        log_scr:str,
        region_file:str=''
    ) -> None:
    """
    Runs bedtools coverage for a list of regions (optional).
    """
    # Make the output directory if it does not exist
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    # Define output file
    output_file = os.path.join(
        output_dir,
        sra_id+'_bedtools_coverage_hist.txt'
    )
    # Initialize the script to run
    l = 'bedtools coverage'
    # Add region file if provided
    if region_file!='':
        l += f' -a {region_file}'
    # Add bam file
    l += f' -b {bam_sorted}'
    # Define as histogram
    l += ' -hist'
    # Define output
    l += f' > {output_file}'
    # Run the script
    log_code(l, log_file=log_scr)
    os.system(l)
    return None

def run_fastqc(
        fastq_file:str,
        output_dir:str,
        log_scr:str
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
    log_code(l, log_file=log_scr)
    os.system(l)
    return None

def run_fastqc_folder(
        fastq_folder:str,
        output_dir:str,
        log_scr:str,
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
    log_code(l, log_file=log_scr)
    os.system(l)
    return None

def run_samtools_depth(
        sra_id:str,
        bam_sorted:str,
        output_dir:str,
        log_scr:str,
        region_file:str=''
    ) -> None:
    """
    Runs samtools depth for a list of regions (optional).
    """
    # Make the output directory if it does not exist
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    # Define output file
    output_file = os.path.join(
        output_dir,
        sra_id+'_samtools_depth.txt'
    )
    # Initialize the script to run
    l = 'samtools depth'
    # Add region file if provided
    if region_file!='':
        l += f' -b {region_file}'
    # Add bam file
    l += f' {bam_sorted}'
    # Define output
    l += f' > {output_file}'
    # Run the script
    log_code(l, log_file=log_scr)
    os.system(l)
    return None
