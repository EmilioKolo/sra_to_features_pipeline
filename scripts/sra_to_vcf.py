#!/usr/bin/env python3


from scripts.snippet import remove_file, run_silent
import json
import logging
import os
import requests
import time


def align_bwa(
        output_folder:str,
        sra_id:str,
        reference_genome:str,
        paired_end:bool,
        log_file:str
    ) -> None:
    """
    Runs an alignment with BWA given an sra_id.
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
    logging.info("Checking SAM file content...")
    l = f'head -n 10 {output_sam}'
    run_silent(l, log_file)
    # Check file size
    l = f'ls -lh {output_sam}'
    run_silent(l, log_file)
    return None

def align_bwa_reads(
        l_reads:list[str],
        output_folder:str,
        out_id:str,
        reference_genome:str,
        paired_end:bool,
        log_file:str
    ) -> None:
    """
    Runs an alignment with BWA given a list of fastq reads.
    """
    # Define if the reads are paired end or single end
    if len(l_reads)>=2:
        paired_end = True
        reads_file_r1 = l_reads[0]
        reads_file_r2 = l_reads[1]
    else:
        paired_end = False
        reads_file_single = l_reads[0]
    # Define output files
    output_prefix = f'{output_folder}/{out_id}'
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
    logging.info("Checking SAM file content...")
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
    # Generate variant file statistics
    logging.info("Variant calling statistics.")
    l = f'bcftools stats {compressed_vcf} > {compressed_vcf}.stats'
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
    # Start a run counter
    i = 0
    while i < 5:
        try:
            run_silent(l, log_file)
            # Exit loop if successful
            break
        except Exception as e:
            logging.info(f"Attempt {i+1} failed: {e}")
            # Check if the files were created and delete them if they exist
            if paired_end:
                if os.path.exists(reads_file_r1):
                    remove_file(reads_file_r1)
                if os.path.exists(reads_file_r2):
                    remove_file(reads_file_r2)
            else:
                if os.path.exists(reads_file_single):
                    remove_file(reads_file_single)
            # Increment the attempt counter
            i += 1
            if i == 5:
                w = "Failed to download FASTQ files after 5 attempts."
                logging.info(w)
                return None
            # Wait before retrying
            time.sleep(60)
    # Visualize file sizes
    t = f'Reads downloaded and extracted to:'
    if paired_end:
        t += f' {reads_file_r1} and {reads_file_r2}'
        logging.info(t)
        # Verify file sizes
        logging.info("Checking file sizes.")
        l = f'du -h {reads_file_r1} {reads_file_r2}'
    else:
        t += f' {reads_file_single}'
        logging.info(t)
        # Verify file sizes
        logging.info("Checking file size.")
        l = f'du -h {reads_file_single}'
    run_silent(l, log_file)
    return None

def get_sra_from_ncbi(sra_accession_id: str) -> dict | None:
    """
    Retrieves SRA (Sequence Read Archive) metadata and download links
    from NCBI via the European Nucleotide Archive (ENA) API.

    Args:
        sra_accession_id (str): The SRA accession ID (e.g., 'SRR000001', 
                                'ERR000001', 'DRR000001').

    Returns:
        dict | None: A dictionary containing SRA metadata and download 
                     links if successful, otherwise None.
    """
    # ENA's API endpoint for searching read runs.
    # We request JSON format and specify the fields we want to retrieve.
    # 'fastq_ftp' and 'sra_ftp' provide direct download links.
    # 'limit=1' ensures we only get one result for a specific accession.
    base_url = "https://www.ebi.ac.uk/ena/portal/api/search"
    f = "run_accession,fastq_ftp,sra_ftp,experiment_accession"
    f += ",sample_accession,study_accession,library_name,library_strategy"
    f += ",library_source,library_selection,instrument_platform"
    f += ",instrument_model,base_count,read_count,scientific_name,tax_id"
    params = {
        "result": "read_run",
        "query": f"run_accession={sra_accession_id}",
        "fields": f,
        "format": "json",
        "limit": 1
    }
    w = f"Attempting to retrieve SRA data for: {sra_accession_id}"
    logging.info(w)
    try:
        # Make the HTTP GET request to the ENA API.
        response = requests.get(base_url, params=params)
        # Raise an HTTPError for bad responses (4xx or 5xx)
        response.raise_for_status()

        # Parse the JSON response.
        data = response.json()

        # The ENA API returns a list of results. For a single accession, 
        # it should be a list with one dictionary, or an empty list 
        # if not found.
        if data:
            sra_info = data[0]
            w = f"Successfully retrieved data for {sra_accession_id}."
            logging.info(w)
            return sra_info
        else:
            w = f"No SRA data found for accession ID: {sra_accession_id}."
            w += ' Please check the ID.'
            logging.info(w)
            return None

    except requests.exceptions.HTTPError as http_err:
        w = f"HTTP error occurred: {http_err}"
        w += f" - Status Code: {response.status_code}"
        logging.info(w)
    except requests.exceptions.ConnectionError as conn_err:
        w = f"Connection error occurred: {conn_err}"
        w += " - Unable to connect to ENA API."
        logging.info(w)
    except requests.exceptions.Timeout as timeout_err:
        w = f"Timeout error occurred: {timeout_err}"
        w += " - Request to ENA API timed out."
        logging.info(w)
    except requests.exceptions.RequestException as req_err:
        w = f"An unexpected error occurred during the request: {req_err}"
        logging.info(w)
    except json.JSONDecodeError as json_err:
        w = f"Error decoding JSON response: {json_err}."
        w += f" Response content: {response.text}"
        logging.info(w)
    except Exception as e:
        logging.info(f"An unexpected error occurred: {e}")

    return None

def sam_to_bam(sam_file:str, bam_file:str, log_file:str) -> None:
    """
    Convert sam file to bam using samtools
    """
    logging.info(f"Converting SAM to BAM: {sam_file} -> {bam_file}...")
    # -b: output BAM
    # -S: input is SAM (optional, but good for clarity)
    l = f'samtools view -bS {sam_file} -o {bam_file}'
    run_silent(l, log_file)

    logging.info(f"SAM to BAM conversion complete: {bam_file}")
    # Check file size
    l = f'ls -lh {bam_file}'
    run_silent(l, log_file)
    # Remove the sam file
    remove_file(sam_file)
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
    logging.info(f"Analyzing variants from {vcf_file} with snpEff...")
    l = f'java -Xmx8g -jar {snpeff_dir}/snpEff/snpEff.jar'
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
    logging.info(f"Sorting BAM file: {bam_file} -> {bam_sorted}...")
    # -o: Output file
    l = f'samtools sort {bam_file} -o {bam_sorted}'
    run_silent(l, log_file)

    logging.info(f"BAM sorting complete: {bam_sorted}")
    # Check file size
    l = f'ls -lh {bam_sorted}'
    run_silent(l, log_file)

    logging.info(f"Indexing sorted BAM file: {bam_sorted}...")
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
    # Remove the bcf file
    remove_file(bcf_file)
    return None
