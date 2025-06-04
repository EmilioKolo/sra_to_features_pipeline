#!/usr/bin/env python3


from collections import defaultdict, Counter
from os import environ
from scripts.feature_generation import extract_regions, fragment_lengths
from scripts.feature_generation import variants_per_bin
from scripts.snippet import run_silent
from scripts.sra_to_vcf import align_bwa, compress_index_vcf
from scripts.sra_to_vcf import download_fastq, get_sra_from_ncbi
from scripts.sra_to_vcf import sam_to_bam, snpeff_analysis
from scripts.sra_to_vcf import sort_index_bam, varcall_mpileup
import logging
import os
import pandas as pd
import sys

# Define logging print level
# Options: INFO, DEBUG, WARNING, ERROR or CRITICAL
logging.basicConfig(level=logging.INFO)

if len(sys.argv) > 1:
    sra_id = sys.argv[1]
    logging.info(f"The provided SRA ID is: {sra_id}")
else:
    logging.info("No SRA ID provided.")
    raise ValueError

# Base values and variables
dict_features = {}
HOME_DIR = os.path.expanduser("~")
output_dir = os.path.join(HOME_DIR, "content/data")
tmp_folder = f'{output_dir}/tmp_{sra_id}'
log_dir = f'{output_dir}/logs'
log_file = f'{log_dir}/{sra_id}'

# Define fasta reference and corresponding gff file
fasta_file = f'{output_dir}/reference.fasta'
gff_file = f'{output_dir}/reference.gff'
reference_genome = fasta_file

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
output_vcf = f"{output_no_tmp}/{sra_id}.vcf"
compressed_vcf = f'{output_vcf}.gz'

# Variables for snpEff
snpeff_dir = f'{output_dir}/bin'
genome_name = 'custom_ref'
snpeff_vcf = f"{output_prefix}.snpeff.vcf"
compressed_snpeff_vcf = f"{snpeff_vcf}.gz"

# Variables for bin creation
genome_sizes = f'{output_dir}/genome.sizes'
bin_size_gvs = 100*1000

# Variables for regions of interest
bed_file = f'regions.bed'
tmp_output = f'{tmp_folder}/counts.csv'

# Variables to be used for Synonymous/Nonsynonymous variant proportion
# vcf file used (must be snpeff, may be compressed)
vcf_file = compressed_snpeff_vcf
# bed files
bed_variants = f'{output_dir}/variants.bed'
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

# Obtain SRA data
sra_data = get_sra_from_ncbi(sra_id)
# Check for files
l_ftp = sra_data['fastq_ftp'].split(';')
if len(l_ftp)==1: # Single end
    paired_end = False
elif len(l_ftp)==2: # Paired end
    paired_end = True
else:
    logging.info(f'WARNING: More than two elements in l_ftp. {l_ftp}')
    paired_end = True

# Download the fastq files from SRA ID
download_fastq(tmp_folder, sra_id, paired_end, log_file)
# Align the reads to the reference genome
align_bwa(tmp_folder, sra_id, reference_genome, paired_end, log_file)
# Convert sam file to bam
sam_to_bam(output_sam, output_bam, log_file)
# Sort and index the bam file
sort_index_bam(output_bam, sorted_bam, log_file)

t = 'Alignment and processing complete.'
t += f' Output BAM: {output_prefix}.sorted.bam'
logging.info(t)
logging.info(f"BAM index: {output_prefix}.sorted.bam.bai")

logging.info("\nAlignment statistics...")
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

# Obtain number of variants in genome bins
bins_dict = variants_per_bin(output_vcf, genome_sizes, bin_size_gvs)
# Load items from bins_dict to dict_features
for key, value in bins_dict.items():
    dict_features[key] = value

# Obtain fragment lengths and their mean, median and st. deviation

# Obtain fragment lengths with samtools (only works with paired ends)
if paired_end:
    fl_mean, fl_median, fl_stdv = fragment_lengths(
        output_dir,
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


# Count variants per region/gene in selected regions/genes

logging.info('Starting to obtain variants per region and genes...')

# Get regions from bed_file
regions = extract_regions(bed_file)

# Count variants per region

# Clear intermediary output file before appending
with open(tmp_output, "w") as f:
    f.write("")

# Load regions and counts into tmp_output
for chrom, start, end, name in regions:
    # Define region string
    region_str = f"{chrom}:{start}-{end}"
    # Add region name
    with open(tmp_output, "a+") as f:
        f.write(f"{name}\t")
    # Use environ to define variables
    environ['region_str'] = region_str
    environ['vcf_file'] = compressed_vcf
    environ['tmp_output'] = tmp_output
    # Use bcftools to check how many variables are in region_str
    l = 'bcftools view -r "$region_str" "$vcf_file" | grep -vc "^#" >> "$tmp_output"'
    run_silent(l, '')

# Read counts back into Python variable
counts = []
with open(tmp_output) as f:
    for line in f:
        name, count = line.strip().split("\t")
        counts.append((name, int(count)))

# Get selected genes from gene_list.txt
selected_genes = open('gene_list.txt', 'r').read().split('\n')

# Save to dict_features
for name, count in counts:
    key = f'{name}_variant_counts'
    dict_features[key] = count
    if name in selected_genes:
        newkey = f'selected_gene_{name}_variant_counts'
        dict_features[newkey] = count


# Synonymous/Nonsynonymous variant proportion per gene

logging.info('Starting to obtain Syn/Nonsyn variant proportion per gene...')

# Define variants using environ
environ['bed_variants'] = bed_variants
environ['vcf_file'] = vcf_file
environ['bed_genes'] = bed_genes
environ['bed_intersect'] = bed_intersect

# Generate bed_variants file
l = 'bcftools query -f \'%CHROM\t%POS\t%END\t%INFO/ANN\n\' $vcf_file | awk \'BEGIN{OFS=\"\t\"} {print $1, $2-1, $2, $4}\' > $bed_variants'
run_silent(l, '')
logging.info(f'{bed_variants} created.')

l = f'bedtools intersect -a {bed_variants} -b {bed_genes} -wa -wb > {bed_intersect}'
run_silent(l, '')
logging.info(f'{bed_intersect} created.')

# Define effect categories
synonymous_terms = {"synonymous_variant"}
nonsynonymous_terms = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "protein_altering_variant",
    "inframe_insertion",
    "inframe_deletion",
    "frameshift_variant"
}

def parse_ann_field(ann_field:str) -> list[str]:
    """
    Parse ANN field and return a list of effects
    """
    # Split the raw annotation field per effect
    annotations = ann_field.split(",")
    effects = []
    # Go through the different variant effects
    for ann in annotations:
        # Split into the different fields
        fields = ann.split("|")
        if len(fields) > 1:
            # Keep the effect field (0 is variant)
            effect = fields[1].strip()
            # Add to output list
            effects.append(effect)
    return effects

# Initialize defaultdict (dictionary) for counts
gene_counts = defaultdict(lambda: {"synonymous": 0, "nonsynonymous": 0})

# Open the intersect bed file
with open(bed_intersect, 'r') as f:
    # Go through lines
    for line in f:
        # Split line
        cols = line.strip().split("\t")
        # Get the gene and annotation fields
        ann_field = cols[3]
        gene = cols[7]
        # Parse the annotationfield
        effects = parse_ann_field(ann_field)
        # Go through obtained effects
        for effect in effects:
            if effect in synonymous_terms:
                gene_counts[gene]["synonymous"] += 1
                break # only count once per variant
            elif effect in nonsynonymous_terms:
                gene_counts[gene]["nonsynonymous"] += 1
                break # only count once per variant

# Log and save results
for gene, counts in gene_counts.items():
    syn = counts["synonymous"]
    nonsyn = counts["nonsynonymous"]
    ratio = syn / nonsyn if nonsyn > 0 else "Inf"
    logging.info(f"{gene}\tSyn: {syn}\tNonsyn: {nonsyn}\tRatio: {ratio}")
    key = f'{gene}_ratio_dN_dS'
    dict_features[key] = ratio


# CNV calling

logging.info('Starting CNV calling...')

# Construct full paths
BAM_PATH = sorted_bam
VCF_PATH = snpeff_vcf
SAMPLE_NAME = sra_id
OUTPUT_DIR = f"{output_prefix}/cnvpytor_results" # Specific output dir for this sample

# Create the output directory if it doesn't exist
l = f'mkdir -p {OUTPUT_DIR}'
run_silent(l, log_file)
logging.info(f"Output directory for results: {OUTPUT_DIR}")

# Define the root file for CNVpytor
ROOT_FILE = os.path.join(OUTPUT_DIR, f"{SAMPLE_NAME}.pytor")
logging.info(f"CNVpytor root file will be: {ROOT_FILE}")

if os.path.exists(ROOT_FILE):
    l = f'rm -r {ROOT_FILE}'
    run_silent(l, '')

# Check if input files exist
if not os.path.exists(BAM_PATH):
    logging.info(f"Error: BAM/CRAM file not found at {BAM_PATH}")
    exit()
if not (os.path.exists(f"{BAM_PATH}.bai") or os.path.exists(f"{BAM_PATH}.crai")):
    logging.info(f"Warning: BAM/CRAM index file not found. Please ensure '{BAM_PATH}.bai' or '{BAM_PATH}.crai' exists alongside the BAM/CRAM. CNVpytor might fail.")

USE_BAF = True # Set to False if you don't have a VCF/GVCF or don't want to use BAF
if USE_BAF and not os.path.exists(VCF_PATH):
    logging.info(f"Warning: VCF/GVCF file not found at {VCF_PATH}. BAF analysis will be skipped.")
    USE_BAF = False
elif USE_BAF and not (os.path.exists(f"{VCF_PATH}.tbi") or os.path.exists(f"{VCF_PATH}.gz.tbi")):
    logging.info(f"Warning: VCF/GVCF index file not found. Ensure '{VCF_PATH}.tbi' or '{VCF_PATH}.gz.tbi' exists alongside the VCF/GVCF. BAF analysis might fail or be inaccurate.")


# CNVpytor Parameters
# Smaller bins offer more resolution, larger bins offer more robustness and detect larger events.
# It's recommended to use a range.
#BIN_SIZES = "1000 10000 100000" # Example: 1kb, 10kb, 100kb bins
BIN_SIZES = str(cnv_bin_size)

# Create root file
l = f'cnvpytor -root {ROOT_FILE} -his {BIN_SIZES} -bam {BAM_PATH}'
run_silent(l, '')

# Process Read Depth (RD) Data
logging.info("\n3. Processing Read Depth (RD) data...")
l = f'cnvpytor -root {ROOT_FILE} -rd {BAM_PATH}'
run_silent(l, '')
logging.info("Read Depth processing complete.")


# Process BAF (BAF) Data
if USE_BAF:
    logging.info("\n4. Processing B-allele Frequency (BAF) data...")
    # First, add SNPs from VCF to the root file
    logging.info(f"Running: cnvpytor -root \"{ROOT_FILE}\" -snp \"{VCF_PATH}\" -sample \"{SAMPLE_NAME}\"")
    l = f'cnvpytor -root {ROOT_FILE} -snp {VCF_PATH} -sample {SAMPLE_NAME}'
    run_silent(l, '')
    # Then, perform BAF analysis with specified bin sizes
    logging.info(f"Running: cnvpytor -root \"{ROOT_FILE}\" -baf {BIN_SIZES}")
    l = f'cnvpytor -root {ROOT_FILE} -baf {BIN_SIZES}'
    run_silent(l, '')
    logging.info("B-allele Frequency processing complete.")
else:
    logging.info("\n4. BAF analysis skipped as specified or due to missing VCF/index.")


# Create Histograms and Partitioning
logging.info("\n5. Generating histograms and partitioning data...")

# Create histograms for RD and BAF (if used)
logging.info(f"Running: cnvpytor -root \"{ROOT_FILE}\" -his {BIN_SIZES}")
l = f'cnvpytor -root {ROOT_FILE} -his {BIN_SIZES} --verbose debug'
run_silent(l, '')

# Partition data for CNV calling
logging.info(f"Running: cnvpytor -root \"{ROOT_FILE}\" -partition {BIN_SIZES}")
l = f'cnvpytor -root {ROOT_FILE} -partition {BIN_SIZES} --verbose debug'
run_silent(l, '')

logging.info("Histograms and partitioning complete.")


# Call CNVs

cnv_call_file = f'{output_prefix}/cnv_calls_{BIN_SIZES}.txt'
logging.info(f"Running: cnvpytor -root \"{ROOT_FILE}\" -call {BIN_SIZES}")
l = f'cnvpytor -root {ROOT_FILE} -call {BIN_SIZES} > {cnv_call_file}'
run_silent(l, log_file)
logging.info("CNV calling complete.")

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

logging.debug(f'{sra_id} features:')
for key, value in dict_features.items():
    logging.debug(f'{key}: {value}')

df_features = pd.DataFrame([dict_features])
try:
    df_features.transpose().to_csv(f'{output_dir}/{sra_id}_features.csv', sep=';', header=[sra_id])
except:
    df_features.transpose().to_csv(f'{output_dir}/{sra_id}_features.csv', sep=';', header=False)

logging.info('Features loaded and saved.\nRemoving temporary files...')

l = f'rm -r {tmp_folder}'
run_silent(l, log_file)

logging.info('Temporary files removed.')
