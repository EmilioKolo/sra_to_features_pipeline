#!/usr/bin/env python3


from collections import defaultdict
from os import environ
from pybedtools import BedTool
import os
import statistics


# Base variables
sra_id = "SRR16993732"
dict_features = {}
output_dir = "/home/emilio/content/data"

# Define fasta reference and corresponding gff file
fasta_file = '/home/emilio/content/data/reference.fasta'
gff_file = '/home/emilio/content/data/reference.gff'
reference_genome = fasta_file

# Variables for fastq files
reads_file_r1 = f"{output_dir}/{sra_id}_1.fastq.gz"
reads_file_r2 = f"{output_dir}/{sra_id}_2.fastq.gz"

# Variables for sam/bam files
output_prefix = f"{output_dir}/{sra_id}"
output_sam = f"{output_prefix}.sam"
output_bam = f"{output_prefix}.bam"
sorted_bam =f"{output_prefix}.sorted.bam"

# Variables for bcf/vcf files
output_bcf = f"{output_prefix}.bcf"
output_vcf = f"{output_prefix}.vcf"
compressed_vcf = f'{output_vcf}.gz'

# Variables for snpEff
snpeff_dir = '/home/emilio/content/data/bin'
genome_name = 'custom_ref'
snpeff_vcf = f"{output_prefix}.snpeff.vcf"
compressed_snpeff_vcf = f"{snpeff_vcf}.gz"

# Variables for bin creation
genome_sizes = "genome.sizes"
bin_size_gvs = 100*1000

# Variables for regions of interest
bed_file = '/home/emilio/content/data/regions.bed'
tmp_output = '/home/emilio/content/data/tmp/counts.csv'

# Variables to be used for Synonymous/Nonsynonymous variant proportion
vcf_file = compressed_snpeff_vcf # vcf file used (must be snpeff, may be compressed)
# bed files
bed_variants = '/home/emilio/content/data/variants.bed'
bed_intersect = '/home/emilio/content/data/intersect.bed'
bed_genes = '/home/emilio/content/data/only_genes.bed'

# Create base output directory
l = f'mkdir -p {output_dir}'
os.system(l)

# Extract fastq files
print(f"\nDownloading and extracting FASTQ files for {sra_id} using fastq-dump...")
# fastq-dump options:
# --split-files: Creates _1.fastq and _2.fastq for paired-end reads
# --gzip: Compresses the output FASTQ files
# -O: Output directory
l = f'fastq-dump --split-files --gzip -O {output_dir} {sra_id}'
os.system(l)
print(f"Reads downloaded and extracted to: {reads_file_r1} and {reads_file_r2}")

# Verify file sizes
print("\nChecking file sizes:")
l = f'du -h {reads_file_r1} {reads_file_r2}'
os.system(l)


# Run alignment
print(f"Aligning reads to {reference_genome} and processing output...")
# -M: Mark shorter split hits as secondary (recommended for Picard compatibility)
# -t: Number of threads (Colab generally has 2 CPU cores available for free tier)
# The '|' pipes the SAM output of bwa to samtools view for conversion to BAM
# samtools sort sorts the BAM file
# samtools index creates the .bai index for quick access
l = f'bwa mem -M -t 2 {reference_genome} {reads_file_r1} {reads_file_r2} > {output_sam}'
os.system(l)

print(f"Alignment to SAM file complete: {output_sam}")

# Check if the SAM file was created and has content
print("\nChecking SAM file content (first 20 lines):")
l = f'head -n 20 {output_sam}'
os.system(l)
# Check file size
l = f'ls -lh {output_sam}'
os.system(l)

print(f"\nConverting SAM to BAM: {output_sam} -> {output_bam}...")
# -b: output BAM
# -S: input is SAM (optional, but good for clarity)
l = f'samtools view -bS {output_sam} -o {output_bam}'
os.system(l)

print(f"SAM to BAM conversion complete: {output_bam}")
# Check file size
l = f'ls -lh {output_bam}'
os.system(l)

print(f"\nSorting BAM file: {output_bam} -> {sorted_bam}...")
# -o: Output file
l = f'samtools sort {output_bam} -o {sorted_bam}'
os.system(l)

print(f"BAM sorting complete: {sorted_bam}")
# Check file size
l = f'ls -lh {sorted_bam}'
os.system(l)

print(f"\nIndexing sorted BAM file: {sorted_bam}...")
l = f'samtools index {sorted_bam}'
os.system(l)

print(f"BAM indexing complete. Index file: {sorted_bam}.bai")
# Check index file size
l = f'ls -lh {sorted_bam}.bai'
os.system(l)

print(f"\nAlignment and processing complete. Output BAM: {output_prefix}.sorted.bam")
print(f"BAM index: {output_prefix}.sorted.bam.bai")

print("\nAlignment statistics:")
l = f'samtools flagstat {output_prefix}.sorted.bam'
os.system(l)


# Run variant calling
print(f"Generating pileup and BCF file for {sorted_bam}...")
# samtools mpileup:
# -u: Uncompressed BCF output (optimal for piping)
# -f: Reference genome file
# -o: Output file
l = f'bcftools mpileup -f {reference_genome} {sorted_bam} > {output_bcf}'
os.system(l)

print(f"Pileup and BCF file generated: {output_bcf}")

# Check bcf file
l = f'tail {output_bcf}'
os.system(l)


# Do calls with BCFtools
print(f"Calling variants from {output_bcf}...")
# bcftools call:
# -m: Multiallelic caller (recommended for most cases)
# -v: Output only variant sites (not homozygous reference sites)
# -o: Output file
l = f'bcftools call -mv -o {output_vcf} {output_bcf}'
os.system(l)

print(f"Variant calling complete. Output VCF: {output_vcf}")
l = f'tail {output_vcf}'
os.system(l)

# Compress and index VCF
print(f"Compressing {output_vcf} with bgzip...")
l = f'bgzip -c {output_vcf} > {compressed_vcf}'
os.system(l)
print(f"Compressed VCF: {compressed_vcf}")

print(f"Indexing {compressed_vcf} with tabix...")
l = f'tabix -p vcf {compressed_vcf}'
os.system(l)
print(f"VCF index: {compressed_vcf}.tbi")

# View the first 50 lines of the VCF (header + some variants)
print("\nFirst 50 lines of the VCF file:")
l = f'zcat {compressed_vcf} | tail'
os.system(l)

print("\nVariant calling statistics:")
l = f'bcftools stats {compressed_vcf} > {output_vcf}.stats'
os.system(l)
l = f'cat {output_vcf}.stats'
os.system(l)


# Analyze variants in VCF with snpeff
l = f'java -Xmx4g -jar {snpeff_dir}/snpEff/snpEff.jar {genome_name} {output_vcf} > {snpeff_vcf}'
os.system(l)
l = f'tail {snpeff_vcf}'
os.system(l)

print(f"Compressing {snpeff_vcf} with bgzip...")
l = f'bgzip -c {snpeff_vcf} > {compressed_snpeff_vcf}'
os.system(l)
print(f"Compressed VCF: {compressed_snpeff_vcf}")

print(f"Indexing {compressed_snpeff_vcf} with tabix...")
l = f'tabix -p vcf {compressed_snpeff_vcf}'
os.system(l)
print(f"VCF index: {compressed_snpeff_vcf}.tbi")


# Obtain number of variants in genome bins

# Generate genome bins
genome = BedTool().window_maker(g=genome_sizes, w=bin_size_gvs)
variants = BedTool(output_vcf)
# Intersect genome and variants to obtain variants per bin
counts = genome.intersect(variants, c=True)
# Put counts in dataframe
l_names = ['chrom', 'start', 'end', 'variant_count']
df_vg_bins = counts.to_dataframe(names=l_names)

# Format dataframe to save as feature
df_vg_bins['bin_region'] = df_vg_bins.apply(lambda row: row['chrom'] + ':' + str(row['start']) + '-' + str(row['end']), axis=1)

total_varcount = sum(df_vg_bins['variant_count'])

for row in df_vg_bins.itertuples():
    key = f'variants_in_{row.bin_region}'
    dict_features[key] = row.variant_count
    key = f'variants_in_{row.bin_region}_normalized'
    dict_features[key] = float(row.variant_count)/total_varcount


# Obtain fragment lengths and their mean, median and st. deviation

# Define environment variables
environ['sorted_bam'] = sorted_bam
environ['sra_id'] = sra_id

# Obtain fragment lengths with samtools
l = 'samtools view -f 0x2 $sorted_bam | awk \'{if ($9>0 && $9<1000) print $9}\' > \"/content/data/fl_\"$sra_id\".txt\"'
os.system(l)

# Open the created file to obtain values
fragment_lengths = open(f"/content/data/fl_{sra_id}.txt", "r").read().splitlines()

# Convert to integers
fragment_lengths = list(map(int, fragment_lengths))

# Get mean, median and standard deviation of fragment lengths (fl)
fl_mean = statistics.mean(fragment_lengths)
fl_median = statistics.median(fragment_lengths)
fl_stdv = statistics.stdev(fragment_lengths)

print('Fragment mean, median and standard deviation:')
print(fl_mean, '/', fl_median, '/', fl_stdv)

dict_features['fragment_lengths_mean'] = fl_mean
dict_features['fragment_lengths_median'] = fl_median
dict_features['fragment_lengths_stdv'] = fl_stdv


# Count variants per region/gene in selected regions/genes

# Get regions from bed_file
regions = []
with open(bed_file) as f:
    # Go through bed file
    for line in f:
        # Discard empty lines and commented lines
        if line.strip() and not line.startswith("#"):
            # Split the line
            fields = line.strip().split('\t')
            # chr, start and end should be in the first 3 fields
            chrom, start, end = fields[:3]
            # Fourth field can be region name (i.e. for genes)
            if len(fields) > 3:
                name = fields[3]
            else:
                name = f"{chrom}:{start}-{end}"
            # Add to region list
            regions.append((chrom, int(start), int(end), name))


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
    os.system(l)

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

# Define variants using os.environ
environ['bed_variants'] = bed_variants
environ['vcf_file'] = vcf_file
environ['bed_genes'] = bed_genes
environ['bed_intersect'] = bed_intersect

# Generate bed_variants file
l = 'bcftools query -f \'%CHROM\t%POS\t%END\t%INFO/ANN\n\' $vcf_file | awk \'BEGIN{OFS=\"\t\"} {print $1, $2-1, $2, $4}\' > $bed_variants'
os.system(l)
print(bed_variants, 'created.')

l = 'bedtools intersect -a $bed_variants -b $bed_genes -wa -wb > $bed_intersect'
os.system(l)
print(bed_intersect, 'created.')

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

# Print and save results
for gene, counts in gene_counts.items():
    syn = counts["synonymous"]
    nonsyn = counts["nonsynonymous"]
    ratio = syn / nonsyn if nonsyn > 0 else "Inf"
    print(f"{gene}\tSyn: {syn}\tNonsyn: {nonsyn}\tRatio: {ratio}")
    key = f'{gene}_ratio_dN_dS'
    dict_features[key] = ratio


# CNV calling

# Construct full paths
BAM_PATH = sorted_bam
VCF_PATH = snpeff_vcf
SAMPLE_NAME = sra_id
OUTPUT_DIR = f"{output_prefix}/cnvpytor_results" # Specific output dir for this sample

# Create the output directory if it doesn't exist
l = f'mkdir -p {OUTPUT_DIR}'
os.system(l)
print(f"Output directory for results: {OUTPUT_DIR}")

# Define the root file for CNVpytor
ROOT_FILE = os.path.join(OUTPUT_DIR, f"{SAMPLE_NAME}.pytor")
print(f"CNVpytor root file will be: {ROOT_FILE}")

# Check if input files exist
if not os.path.exists(BAM_PATH):
    print(f"Error: BAM/CRAM file not found at {BAM_PATH}")
    exit()
if not (os.path.exists(f"{BAM_PATH}.bai") or os.path.exists(f"{BAM_PATH}.crai")):
    print(f"Warning: BAM/CRAM index file not found. Please ensure '{BAM_PATH}.bai' or '{BAM_PATH}.crai' exists alongside the BAM/CRAM. CNVpytor might fail.")

USE_BAF = True # Set to False if you don't have a VCF/GVCF or don't want to use BAF
if USE_BAF and not os.path.exists(VCF_PATH):
    print(f"Warning: VCF/GVCF file not found at {VCF_PATH}. BAF analysis will be skipped.")
    USE_BAF = False
elif USE_BAF and not (os.path.exists(f"{VCF_PATH}.tbi") or os.path.exists(f"{VCF_PATH}.gz.tbi")):
    print(f"Warning: VCF/GVCF index file not found. Ensure '{VCF_PATH}.tbi' or '{VCF_PATH}.gz.tbi' exists alongside the VCF/GVCF. BAF analysis might fail or be inaccurate.")

# --- CNVpytor Parameters ---
# Smaller bins offer more resolution, larger bins offer more robustness and detect larger events.
# It's recommended to use a range.
BIN_SIZES = "1000 10000 100000" # Example: 1kb, 10kb, 100kb bins
#BIN_SIZES = "1000"

# --- 3. Process Read Depth (RD) Data ---
print("\n3. Processing Read Depth (RD) data...")
l = f'cnvpytor -root {ROOT_FILE} -rd {BAM_PATH}'
os.system(l)
print("Read Depth processing complete.")

# --- 4. Process BAF (BAF) Data ---

if USE_BAF:
    print("\n4. Processing B-allele Frequency (BAF) data...")
    # First, add SNPs from VCF to the root file
    print(f"Running: cnvpytor -root \"{ROOT_FILE}\" -snp \"{VCF_PATH}\" -sample \"{SAMPLE_NAME}\"")
    l = f'cnvpytor -root {ROOT_FILE} -snp {VCF_PATH} -sample {SAMPLE_NAME}'
    os.system(l)
    # Then, perform BAF analysis with specified bin sizes
    print(f"Running: cnvpytor -root \"{ROOT_FILE}\" -baf {BIN_SIZES}")
    l = f'cnvpytor -root {ROOT_FILE} -baf {BIN_SIZES}'
    os.system(l)
    print("B-allele Frequency processing complete.")
else:
    print("\n4. BAF analysis skipped as specified or due to missing VCF/index.")

# --- 5. Create Histograms and Partitioning ---
print("\n5. Generating histograms and partitioning data...")

# Create histograms for RD and BAF (if used)
print(f"Running: cnvpytor -root \"{ROOT_FILE}\" -his {BIN_SIZES}")
l = f'cnvpytor -root {ROOT_FILE} -his {BIN_SIZES} --verbose debug'
os.system(l)

# Partition data for CNV calling
print(f"Running: cnvpytor -root \"{ROOT_FILE}\" -partition {BIN_SIZES}")
l = f'cnvpytor -root {ROOT_FILE} -partition {BIN_SIZES} --verbose debug'
os.system(l)

print("Histograms and partitioning complete.")

# --- 6. Call CNVs ---

print(f"Running: cnvpytor -root \"{ROOT_FILE}\" -call {BIN_SIZES}")
l = f'cnvpytor -root {ROOT_FILE} -call {BIN_SIZES}'
os.system(l)
print("CNV calling complete.")
