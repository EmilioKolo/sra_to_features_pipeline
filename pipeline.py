#!/usr/bin/env python3


from collections import Counter
from scripts.checks import classify_unaligned_reads
from scripts.feature_generation import cnvpytor_pipeline
from scripts.feature_generation import count_syn_nonsyn, create_counts
from scripts.feature_generation import extract_regions, fragment_lengths
from scripts.feature_generation import parse_gff_for_genes
from scripts.feature_generation import variants_per_bin_os
from scripts.snippet import change_output_ownership, define_bname
from scripts.snippet import get_value, run_silent
from scripts.sra_to_vcf import align_bwa, align_bwa_reads
from scripts.sra_to_vcf import compress_index_vcf
from scripts.sra_to_vcf import download_fastq, get_sra_from_ncbi
from scripts.sra_to_vcf import sam_to_bam, snpeff_analysis
from scripts.sra_to_vcf import sort_index_bam, varcall_mpileup
import logging
import os
import pandas as pd
import random
import sys
import time

# Define logging print level
# Options: DEBUG, INFO, WARNING, ERROR or CRITICAL
logging.basicConfig(level=logging.INFO)

if len(sys.argv) > 1:
    if sys.argv[1].endswith(('.fastq', '.fastq.gz')):
        fastq_given = True
        in_path = '/content/data/input'
        if len(sys.argv) > 2:
            reads_file_r1 = os.path.join(in_path, sys.argv[1])
            reads_file_r2 = os.path.join(in_path, sys.argv[2])
            paired_end = True
            # Define basename
            sra_id = define_bname(
                sys.argv[1],
                sys.argv[2]
            ).replace('.', '_')
            w = "The provided fastq files are:"
            w += f' {reads_file_r1} and {reads_file_r2}'
            logging.info(w)
        else:
            reads_file_single = os.path.join(in_path, sys.argv[1])
            paired_end = False
            # Define basename
            sra_id = sys.argv[1].split('.')[0]
            w = f"The provided fastq file is: {reads_file_single}"
            logging.info(w)
        output_dir = '/content/data/output'
        logging.info(f'The output directory is "{output_dir}"')
    else:
        fastq_given = False
        sra_id = sys.argv[1]
        logging.info(f"The provided SRA ID is: {sra_id}")
        if len(sys.argv) > 2:
            output_dir = sys.argv[2].rstrip('/')
            # Use the absolute path for output_dir
            output_dir = os.path.abspath(output_dir)
            logging.info(f'The output directory is "{output_dir}"')
        else:
            output_dir = '/content/data/output'
            w = f'Output directory not provided, using "{output_dir}"'
            logging.info(w)
else:
    logging.info("No SRA ID nor fastq files provided.")
    raise ValueError

print(f'RUNNING FOR {sra_id}')

# Define base directory (for non-output files)
BASE_DIR = get_value('BASE_DIR').strip('"').rstrip('/')
# Base values and variables
dict_features = {}
tmp_folder = f'{BASE_DIR}/tmp_{sra_id}'
log_dir = f'{output_dir}/logs'
log_file = f'{log_dir}/{sra_id}'
THREADS = int(get_value('THREADS'))

# Define fasta reference and corresponding gff file
fasta_file = f'{BASE_DIR}/data/reference.fasta'
gff_file = f'{BASE_DIR}/data/reference.gff'
reference_genome = fasta_file

if not fastq_given:
    # Variables for fastq files
    reads_file_r1 = f"{tmp_folder}/{sra_id}_1.fastq.gz"
    reads_file_r2 = f"{tmp_folder}/{sra_id}_2.fastq.gz"
    reads_file_single = f"{tmp_folder}/{sra_id}.fastq.gz"

# Variables for sam/bam files
output_no_tmp = f'{output_dir}/{sra_id}'
output_prefix = f"{tmp_folder}/{sra_id}"
output_sam = f"{output_prefix}.sam"
output_bam = f"{output_prefix}.bam"
sorted_bam =f"{output_prefix}.sorted.bam"

# Variables for bcf/vcf files
output_bcf = f"{output_prefix}.bcf"
output_vcf = f"{output_prefix}.vcf"
compressed_vcf = f'{output_no_tmp}/{sra_id}.vcf.gz'

# Variables for snpEff
snpeff_dir = f'{BASE_DIR}/data/bin'
genome_name = get_value('GENOME_NAME').strip('"')
snpeff_vcf = f"{output_prefix}.snpeff.vcf"
compressed_snpeff_vcf = f"{snpeff_vcf}.gz"

# Variables for kraken2
kraken_db = f'{BASE_DIR}/install/kraken2-db'

# Variables for bin creation
genome_sizes = f'{BASE_DIR}/data/genome.sizes'
bin_size_gvs = 100*1000

# Variables for regions of interest
bed_file = f'regions.bed'
counts_file = f'{tmp_folder}/counts.csv'

# Variables to be used for Synonymous/Nonsynonymous variant proportion
# vcf file used (must be snpeff, may be compressed)
vcf_file = compressed_snpeff_vcf
# bed files
bed_variants = f'{tmp_folder}/variants.bed'
bed_intersect = f'{tmp_folder}/intersect.bed'
bed_genes = f'{output_dir}/only_genes.bed'

# Bin sizes for CNVpytor (recommended: 1k, 10k, 100k)
cnv_bin_size = 100*1000

# Create base output directories
l = f'mkdir -p {output_dir}'
run_silent(l, '')
l = f'mkdir -p {tmp_folder}'
run_silent(l, '')
l = f'mkdir -p {log_dir}'
run_silent(l, '')
l = f'mkdir -p {output_no_tmp}'
run_silent(l, '')

if not fastq_given:
    # Obtain SRA data
    sra_data = get_sra_from_ncbi(sra_id)
    cont = 0
    while sra_data is None:
        # Add random wait time
        time.sleep(random.randint(15, 30))
        w = f"Retrying to get SRA data for {sra_id} (attempt {cont+1})"
        logging.info(w)
        # Increment the counter
        cont+=1
        # Retry getting SRA data
        sra_data = get_sra_from_ncbi(sra_id)
        if cont>5:
            er = f"Could not retrieve SRA data for ID: {sra_id}"
            logging.error(er)
            sys.exit(1)
    # Check for files
    l_ftp = sra_data['fastq_ftp'].split(';')
    if len(l_ftp)==1: # Single end
        paired_end = False
    elif len(l_ftp)==2: # Paired end
        paired_end = True
    else:
        logging.info(f'WARNING: More than two elements in l_ftp.\n{l_ftp}')
        paired_end = True

######################### SRA TO VCF PIPELINE ###########################

if not fastq_given:
    # Download the fastq files from SRA ID
    download_fastq(tmp_folder, sra_id, paired_end, log_file)
# Align the reads to the reference genome
if fastq_given:
    if paired_end:
        l_reads = [reads_file_r1, reads_file_r2]
    else:
        l_reads = [reads_file_single]
    align_bwa_reads(
        l_reads,
        tmp_folder,
        sra_id,
        reference_genome,
        paired_end,
        log_file
    )
else:
    align_bwa(tmp_folder, sra_id, reference_genome, paired_end, log_file)
# Convert sam file to bam
sam_to_bam(output_sam, output_bam, log_file)
# Sort and index the bam file
sort_index_bam(output_bam, sorted_bam, log_file)
# Finishing messages
t = 'Alignment and processing complete.'
t += f' Output BAM: {output_prefix}.sorted.bam'
logging.info(t)
logging.info(f"BAM index: {output_prefix}.sorted.bam.bai")
logging.info("Alignment statistics...")
l = f'samtools flagstat {output_prefix}.sorted.bam'
run_silent(l, log_file)

# Run variant calling
varcall_mpileup(
    sorted_bam,
    output_bcf,
    output_vcf,
    reference_genome,
    log_file
)
# Compress and index the created vcf file
compress_index_vcf(output_vcf, compressed_vcf, log_file)
# Analyze variants in VCF with snpeff
snpeff_analysis(
    output_vcf,
    snpeff_vcf,
    compressed_snpeff_vcf,
    genome_name,
    snpeff_dir,
    log_file
)

################################ CHECKS #################################

# Check unaligned reads with Kraken2
try:
    logging.info("Classifying unaligned reads with Kraken2...")
    base_kraken_out = f'{output_no_tmp}/{sra_id}'
    classify_unaligned_reads(
        output_bam,
        base_kraken_out,
        kraken_db=kraken_db,
        threads=THREADS
    )
    logging.info("Unaligned reads classified successfully.")
    p = f"Results in: {output_no_tmp}_kraken2_output.txt"
    p += f' and {output_no_tmp}_kraken2_report.txt'
    logging.info(p)
except Exception as e:
    logging.error(f"Error classifying unaligned reads: {e}")
    logging.info("Skipping unaligned reads classification.")

########################## FEATURE GENERATION ###########################

# Obtain fragment lengths and their mean, median and st. deviation
if paired_end:
    # Obtain fragment lengths with samtools (only works with paired ends)
    fl_mean, fl_median, fl_stdv = fragment_lengths(
        tmp_folder,
        sra_id,
        sorted_bam
    )
else:
    fl_mean = 'NA'
    fl_median = 'NA'
    fl_stdv = 'NA'
    logging.info('Single-end reads. Fragment length not calculated.')
# Assign values to dict_features
dict_features['fragment_lengths_mean'] = fl_mean
dict_features['fragment_lengths_median'] = fl_median
dict_features['fragment_lengths_stdv'] = fl_stdv

# Obtain number of variants in genome bins
bins_dict = variants_per_bin_os(output_vcf, genome_sizes, bin_size_gvs)
# Load items from bins_dict to dict_features
for key, value in bins_dict.items():
    dict_features[key] = value

# Count variants per region/gene in selected regions/genes
logging.info('Starting to obtain variants per region and genes...')
# Get regions from bed_file
selected_regions = extract_regions(bed_file)
all_genes = parse_gff_for_genes(gff_file)
regions = selected_regions + all_genes
# Create counts file
create_counts(
    counts_file,
    regions,
    compressed_vcf
)
# Read counts back into Python variable
counts = []
with open(counts_file) as f:
    for line in f:
        name, count = line.strip().split("\t")
        counts.append((name, int(count)))
# Get selected genes from gene_list.txt
selected_genes = open('gene_list.txt', 'r').read().split('\n')
# Save to dict_features
for name, count in counts:
    if name in selected_genes:
        newkey = f'selected_gene_{name}_variant_counts'
        dict_features[newkey] = count
    else:
        key = f'{name}_variant_counts'
        dict_features[key] = count

# Synonymous/Nonsynonymous variant proportion per gene
gene_counts = count_syn_nonsyn(
    bed_variants,
    vcf_file,
    bed_genes,
    bed_intersect
)
# Log and save results
for gene, counts in gene_counts.items():
    syn = counts["synonymous"]
    nonsyn = counts["nonsynonymous"]
    ratio = syn / nonsyn if nonsyn > 0 else "Inf"
    w = f"{gene}\tSyn: {syn}\tNonsyn: {nonsyn}\tRatio: {ratio}"
    logging.info(w)
    key = f'{gene}_ratio_dN_dS'
    dict_features[key] = ratio

# CNV calling
cnv_call_file = cnvpytor_pipeline(
    output_prefix,
    sra_id,
    sorted_bam,
    snpeff_vcf,
    cnv_bin_size,
    log_file
)
# Read CNVpytor output with python 
calls = []
with open(cnv_call_file, 'r') as f:
    for line in f.readlines():
        calls.append(line.rstrip('\n').split('\t'))
# Count CNVs per chromosome
cnv_counts = Counter()
for call in calls:
    chrom = call[1].split(':')[0]
    cnv_counts[chrom] += 1
# Save result
for chrom, count in sorted(cnv_counts.items()):
    if str(chrom).startswith('chr'):
        logging.info(f"{chrom}: {count} CNVs")
        key = f'{chrom}_cnv_count'
    else:
        logging.info(f"chr{chrom}: {count} CNVs")
        key = f'chr{chrom}_cnv_count'
    dict_features[key] = count

################### SAVING FEATURES / DELETING FILES ####################

df_features = pd.DataFrame([dict_features])
df_out = f'{output_dir}/{sra_id}_features.csv'
try:
    df_features.transpose().to_csv(df_out, sep=';', header=[sra_id])
except:
    df_features.transpose().to_csv(df_out, sep=';', header=False)

logging.info('Features loaded and saved.')
logging.info('Removing temporary files...')

# Remove temporaryd folder
l = f'rm -r {tmp_folder}'
run_silent(l, log_file)

w = 'Temporary files removed. Changing ownership of output files...'
logging.info(w)
# Change ownership of output files to the host user
change_output_ownership(output_dir)
