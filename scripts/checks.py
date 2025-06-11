#!/usr/bin/env python3


import logging
import os


def classify_unaligned_reads(
        bam_file:str,
        output_prefix:str,
        kraken_db:str="~/kraken2-db",
        threads:int=4
        ) -> None:
    """
    Classify unaligned reads from a BAM file using Kraken2.
    """
    # Define created file paths
    unaligned_bam = f"{output_prefix}_unaligned.bam"
    unaligned_fastq = f"{output_prefix}_unaligned.fastq"
    kraken_output = f"{output_prefix}_kraken2_output.txt"
    kraken_report = f"{output_prefix}_kraken2_report.txt"
    # Extract unaligned reads
    logging.info('Extracting unaligned reads from BAM file...')
    l = f'samtools view -f 4 -b {bam_file} -o {unaligned_bam}'
    os.system(l)
    logging.info('Converting unaligned BAM to FASTQ...')
    l = f'samtools fastq {unaligned_bam} -o {unaligned_fastq}'
    os.system(l)
    logging.info('Running Kraken2 classification...')
    # Classify with Kraken2
    l = f'kraken2 --db {kraken_db} --threads {threads}'
    l += f' --output {kraken_output} --report {kraken_report}'
    l += f' {unaligned_fastq}'
    os.system(l)
    p = 'Kraken2 classification complete. Results in:'
    p += f'\n * {kraken_output}\n * {kraken_report}'
    print(p)
    # Clean up intermediate files
    os.remove(unaligned_bam)
    os.remove(unaligned_fastq)
    return None
