#!/usr/bin/env python3


from scripts.snippet import run_silent
import json
import logging
import requests


def align_bwa(
        output_folder:str,
        sra_id:str,
        reference_genome:str,
        paired_end:bool,
        log_file:str
        ) -> None:
    """
    Runs an alignment with BWA.
    """
    # Define file names
    reads_file_r1 = f'{output_folder}/{sra_id}_1.fastq.gz'
    reads_file_r2 = f'{output_folder}/{sra_id}_2.fastq.gz'
    reads_file_single = f'{output_folder}/{sra_id}.fastq.gz'
    output_prefix = f'{output_folder}/{sra_id}'
    output_sam = f'{output_prefix}.sam'
    # Run alignment
    t = f'Aligning reads to {reference_genome}'
    t += ' and processing output...'
    logging.info(t)
    # -M: Mark shorter split hits as secondary (recommended for Picard compatibility)
    # -t: Number of threads (Colab generally has 2 CPU cores available for free tier)
    # The '|' pipes the SAM output of bwa to samtools view for conversion to BAM
    # samtools sort sorts the BAM file
    # samtools index creates the .bai index for quick access
    l = 'bwa mem -M -t 2'
    l += f' {reference_genome}'
    if paired_end:
        l += f' {reads_file_r1} {reads_file_r2}'
    else:
        l += f' {reads_file_single}'
    l += f' > {output_sam}'
    run_silent(l, log_file)
    logging.info(f"Alignment to SAM file complete: {output_sam}")
    # Check if the SAM file was created and has content
    logging.info("\nChecking SAM file content...")
    l = f'head -n 10 {output_sam}'
    run_silent(l, log_file)
    # Check file size
    l = f'ls -lh {output_sam}'
    run_silent(l, log_file)
    return None

def compress_index_vcf(
        vcf_file:str,
        compressed_vcf:str,
        log_file:str
        ) -> None:
    """
    Compresses a vcf file, indexes it with tabix and prints some stats.
    """
    # Compress VCF
    logging.info(f"Compressing {vcf_file} with bgzip...")
    l = f'bgzip -c {vcf_file} > {compressed_vcf}'
    run_silent(l, log_file)
    logging.info(f"Compressed VCF: {compressed_vcf}")

    # Index VCF with tabix
    logging.info(f"Indexing {compressed_vcf} with tabix...")
    l = f'tabix -p vcf {compressed_vcf}'
    run_silent(l, log_file)
    logging.info(f"VCF index: {compressed_vcf}.tbi")

    # View the first lines of the VCF (header + some variants)
    logging.info("First lines of the VCF file.")
    l = f'zcat {compressed_vcf} | tail'
    run_silent(l, log_file)

    logging.info("Variant calling statistics.")
    l = f'bcftools stats {compressed_vcf} > {vcf_file}.stats'
    run_silent(l, log_file)
    l = f'cat {vcf_file}.stats'
    run_silent(l, log_file)
    return None

def download_fastq(
        output_folder:str,
        sra_id:str,
        paired_end:bool,
        log_file:str
        ) -> None:
    """
    Download fastq files from an SRA ID using fastq-dump.
    """
    # Define file names
    reads_file_r1 = f'{output_folder}/{sra_id}_1.fastq.gz'
    reads_file_r2 = f'{output_folder}/{sra_id}_2.fastq.gz'
    reads_file_single = f'{output_folder}/{sra_id}.fastq.gz'
    t = 'Downloading and extracting FASTQ files for'
    t += f' {sra_id} using fastq-dump...'
    logging.info(t)
    # fastq-dump options:
    # --gzip: Compresses the output FASTQ files
    # -O: Output directory
    l = 'fastq-dump --gzip'
    if paired_end:
        l += ' --split-files'
    l += f' -O {output_folder} {sra_id}'
    run_silent(l, log_file)
    # Visualize file sizes
    t = f'Reads downloaded and extracted to:'
    if paired_end:
        t += f' {reads_file_r1} and {reads_file_r2}'
        logging.info(t)
        # Verify file sizes
        logging.info("\nChecking file sizes.")
        l = f'du -h {reads_file_r1} {reads_file_r2}'
    else:
        t += f' {reads_file_single}'
        logging.info(t)
        # Verify file sizes
        logging.info("\nChecking file size.")
        l = f'du -h {reads_file_single}'
    run_silent(l, log_file)
    return None

def get_sra_from_ncbi(sra_accession_id: str) -> dict | None:
    """
    Retrieves SRA (Sequence Read Archive) metadata and download links
    from NCBI via the European Nucleotide Archive (ENA) API.

    Args:
        sra_accession_id (str): The SRA accession ID (e.g., 'SRR000001', 'ERR000001', 'DRR000001').

    Returns:
        dict | None: A dictionary containing SRA metadata and download links if successful,
                     otherwise None.
    """
    # ENA's API endpoint for searching read runs.
    # We request JSON format and specify the fields we want to retrieve.
    # 'fastq_ftp' and 'sra_ftp' provide direct download links.
    # 'limit=1' ensures we only get one result for a specific accession.
    base_url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "read_run",
        "query": f"run_accession={sra_accession_id}",
        "fields": "run_accession,fastq_ftp,sra_ftp,experiment_accession,sample_accession,study_accession,library_name,library_strategy,library_source,library_selection,instrument_platform,instrument_model,base_count,read_count,scientific_name,tax_id",
        "format": "json",
        "limit": 1
    }

    logging.info(f"Attempting to retrieve SRA data for: {sra_accession_id}")
    try:
        # Make the HTTP GET request to the ENA API.
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)

        # Parse the JSON response.
        data = response.json()

        # The ENA API returns a list of results. For a single accession, it should be a list
        # with one dictionary, or an empty list if not found.
        if data:
            sra_info = data[0]
            logging.info(f"Successfully retrieved data for {sra_accession_id}.")
            return sra_info
        else:
            logging.info(f"No SRA data found for accession ID: {sra_accession_id}. Please check the ID.")
            return None

    except requests.exceptions.HTTPError as http_err:
        logging.info(f"HTTP error occurred: {http_err} - Status Code: {response.status_code}")
    except requests.exceptions.ConnectionError as conn_err:
        logging.info(f"Connection error occurred: {conn_err} - Unable to connect to ENA API.")
    except requests.exceptions.Timeout as timeout_err:
        logging.info(f"Timeout error occurred: {timeout_err} - Request to ENA API timed out.")
    except requests.exceptions.RequestException as req_err:
        logging.info(f"An unexpected error occurred during the request: {req_err}")
    except json.JSONDecodeError as json_err:
        logging.info(f"Error decoding JSON response: {json_err}. Response content: {response.text}")
    except Exception as e:
        logging.info(f"An unexpected error occurred: {e}")

    return None

def sam_to_bam(sam_file:str, bam_file:str, log_file:str) -> None:
    """
    Convert sam file to bam using samtools
    """
    logging.info(f"\nConverting SAM to BAM: {sam_file} -> {bam_file}...")
    # -b: output BAM
    # -S: input is SAM (optional, but good for clarity)
    l = f'samtools view -bS {sam_file} -o {bam_file}'
    run_silent(l, log_file)

    logging.info(f"SAM to BAM conversion complete: {bam_file}")
    # Check file size
    l = f'ls -lh {bam_file}'
    run_silent(l, log_file)
    return None

def snpeff_analysis(
        vcf_file:str,
        snpeff_vcf:str,
        compressed_snpeff_vcf:str,
        genome_name:str,
        snpeff_dir:str,
        log_file:str
        ) -> None:
    """
    Runs snpeff on a vcf file, then compresses and indexes the output.
    """
    # Analyze variants in VCF with snpeff
    l = f'java -Xmx4g -jar {snpeff_dir}/snpEff/snpEff.jar'
    l += f' {genome_name} {vcf_file} > {snpeff_vcf}'
    run_silent(l, log_file)
    # Visualize the snpeff vcf file
    l = f'tail {snpeff_vcf}'
    run_silent(l, log_file)
    # Run compression with bgzip
    logging.info(f"Compressing {snpeff_vcf} with bgzip...")
    l = f'bgzip -c {snpeff_vcf} > {compressed_snpeff_vcf}'
    run_silent(l, log_file)
    logging.info(f"Compressed VCF: {compressed_snpeff_vcf}")
    # Index the compressed snpeff vcf file
    logging.info(f"Indexing {compressed_snpeff_vcf} with tabix...")
    l = f'tabix -p vcf {compressed_snpeff_vcf}'
    run_silent(l, log_file)
    logging.info(f"VCF index: {compressed_snpeff_vcf}.tbi")
    return None

def sort_index_bam(bam_file:str, bam_sorted:str, log_file:str) -> None:
    """
    Runs samtools sort and samtools index on a bam file.
    """
    logging.info(f"\nSorting BAM file: {bam_file} -> {bam_sorted}...")
    # -o: Output file
    l = f'samtools sort {bam_file} -o {bam_sorted}'
    run_silent(l, log_file)

    logging.info(f"BAM sorting complete: {bam_sorted}")
    # Check file size
    l = f'ls -lh {bam_sorted}'
    run_silent(l, log_file)

    logging.info(f"\nIndexing sorted BAM file: {bam_sorted}...")
    l = f'samtools index {bam_sorted}'
    run_silent(l, log_file)

    logging.info(f"BAM indexing complete. Index file: {bam_sorted}.bai")
    # Check index file size
    l = f'ls -lh {bam_sorted}.bai'
    run_silent(l, log_file)
    return None

def varcall_mpileup(
        bam_sorted:str,
        bcf_file:str,
        vcf_file:str,
        reference_genome:str,
        log_file:str
        ) -> None:
    """
    Perform variant calling using bcftools mpileup.
    """
    # Run variant calling
    logging.info(f"Generating pileup and BCF file for {bam_sorted}...")
    # samtools mpileup:
    # -u: Uncompressed BCF output (optimal for piping)
    # -f: Reference genome file
    # -o: Output file
    l = f'bcftools mpileup -f {reference_genome} {bam_sorted} > {bcf_file}'
    run_silent(l, log_file)

    logging.info(f"Pileup and BCF file generated: {bcf_file}")
    # Check bcf file
    l = f'tail {bcf_file}'
    run_silent(l, log_file)

    # Do calls with BCFtools
    logging.info(f"Calling variants from {bcf_file}...")
    # bcftools call:
    # -m: Multiallelic caller (recommended for most cases)
    # -v: Output only variant sites (not homozygous reference sites)
    # -o: Output file
    l = f'bcftools call -mv -o {vcf_file} {bcf_file}'
    run_silent(l, log_file)

    logging.info(f"Variant calling complete. Output VCF: {vcf_file}")
    l = f'tail {vcf_file}'
    run_silent(l, log_file)
    return None
