#!/usr/bin/env python3


from collections import defaultdict, Counter
from os import environ
from pybedtools import BedTool
from scripts.snippet import run_silent
import logging
import os
import statistics


def cnvpytor_pipeline(
        output_prefix:str,
        sra_id:str,
        bam_sorted:str,
        vcf_snpeff:str,
        bin_size:int,
        log_file:str,
        use_baf:bool = True
        ) -> str:
    logging.info('Starting CNV calling...')
    # Construct full paths
    BAM_PATH = bam_sorted
    VCF_PATH = vcf_snpeff
    SAMPLE_NAME = sra_id
    # Define a specific output dir for this sample
    OUTPUT_DIR = f"{output_prefix}/cnvpytor_results"
    # Create the output directory if it doesn't exist
    l = f'mkdir -p {OUTPUT_DIR}'
    run_silent(l, log_file)
    logging.info(f"Output directory for results: {OUTPUT_DIR}")
    # Define the root file for CNVpytor
    ROOT_FILE = os.path.join(OUTPUT_DIR, f"{SAMPLE_NAME}.pytor")
    logging.info(f"CNVpytor root file will be: {ROOT_FILE}")
    # Remove existing root file if it exists
    if os.path.exists(ROOT_FILE):
        l = f'rm -r {ROOT_FILE}'
        run_silent(l, '')
    # Check if input files exist
    bool_bam = os.path.exists(BAM_PATH)
    bool_bai = os.path.exists(f"{BAM_PATH}.bai")
    bool_crai = os.path.exists(f"{BAM_PATH}.crai")
    if not bool_bam:
        logging.info(f"Error: BAM/CRAM file not found at {BAM_PATH}")
        exit()
    if not (bool_bai or bool_crai):
        t = "Warning: BAM/CRAM index file not found."
        t += f" Please ensure '{BAM_PATH}.bai' or '{BAM_PATH}.crai'"
        t += " exists alongside the BAM/CRAM. CNVpytor might fail."
        logging.info(t)
    # Define if you have a VCF/GVCF and want to use BAF
    USE_BAF = use_baf
    bool_vcf = os.path.exists(VCF_PATH)
    bool_tbi = os.path.exists(f"{VCF_PATH}.tbi")
    bool_gz_tbi = os.path.exists(f"{VCF_PATH}.gz.tbi")
    if USE_BAF and not bool_vcf:
        t = f"Warning: VCF/GVCF file not found at {VCF_PATH}."
        t += " BAF analysis will be skipped."
        logging.info(t)
        USE_BAF = False
    elif USE_BAF and not (bool_tbi or bool_gz_tbi):
        t = "Warning: VCF/GVCF index file not found."
        t += f" Ensure '{VCF_PATH}.tbi' or '{VCF_PATH}.gz.tbi'"
        t += " exists alongside the VCF/GVCF."
        t += " BAF analysis might fail or be inaccurate."
        logging.info(t)
    # CNVpytor Parameters
    # Smaller bins offer more resolution
    # Larger bins offer more robustness and detect larger events
    # It's recommended to use a range.
    #BIN_SIZES = "1000 10000 100000" # Example: 1kb, 10kb, 100kb bins
    BIN_SIZES = str(bin_size)

    # Create root file
    l = f'cnvpytor -root {ROOT_FILE} -his {BIN_SIZES} -bam {BAM_PATH}'
    run_silent(l, '')

    # Process Read Depth (RD) Data
    logging.info("3. Processing Read Depth (RD) data...")
    l = f'cnvpytor -root {ROOT_FILE} -rd {BAM_PATH}'
    run_silent(l, '')
    logging.info("Read Depth processing complete.")

    # Process BAF (BAF) Data
    if USE_BAF:
        logging.info("4. Processing B-allele Frequency (BAF) data...")
        # First, add SNPs from VCF to the root file
        l = f'cnvpytor -root {ROOT_FILE}'
        l += f' -snp {VCF_PATH} -sample {SAMPLE_NAME}'
        run_silent(l, '')
        # Then, perform BAF analysis with specified bin sizes
        l = f'cnvpytor -root {ROOT_FILE} -baf {BIN_SIZES}'
        run_silent(l, '')
        t = "B-allele Frequency processing complete."
        logging.info(t)
    else:
        t = "4. BAF analysis skipped as specified"
        t += " or due to missing VCF/index."
        logging.info(t)
    # Create Histograms and Partitioning
    t = "5. Generating histograms and partitioning data..."
    logging.info(t)
    # Create histograms for RD and BAF (if used)
    l = f'cnvpytor -root {ROOT_FILE} -his {BIN_SIZES} --verbose debug'
    run_silent(l, '')
    # Partition data for CNV calling
    l = f'cnvpytor -root {ROOT_FILE}'
    l += f' -partition {BIN_SIZES} --verbose debug'
    run_silent(l, '')
    logging.info("Histograms and partitioning complete.")
    # Call CNVs
    cnv_call_file = f'{output_prefix}/cnv_calls_{BIN_SIZES}.txt'
    l = f'cnvpytor -root {ROOT_FILE} -call {BIN_SIZES} > {cnv_call_file}'
    run_silent(l, log_file)
    logging.info("CNV calling complete.")
    return cnv_call_file

def create_counts(
        counts_file:str,
        regions:list[list[str,int,int,str]],
        vcf_compressed:str
        ) -> None:
    # Clear intermediary output file before appending
    with open(counts_file, "w") as f:
        f.write("")

    # Load regions and counts into out_folder
    for chrom, start, end, name in regions:
        # Define region string
        region_str = f"{chrom}:{start}-{end}"
        # Add region name
        with open(counts_file, "a+") as f:
            f.write(f"{name}\t")
        # Use environ to define variables
        environ['region_str'] = region_str
        environ['vcf_file'] = vcf_compressed
        environ['tmp_output'] = counts_file
        # Use bcftools to check how many variables are in region_str
        l = 'bcftools view -r "$region_str"'
        l += ' "$vcf_file" | grep -vc "^#" >> "$tmp_output"'
        run_silent(l, '')
    return None

def count_syn_nonsyn(
        bed_variants:str,
        vcf_file:str,
        bed_genes:str,
        bed_intersect:str
        ) -> dict[str, dict[str, int]]:
    """
    Counts the proportion of synonymous to nonsynonymous variants per gene.
    """
    t = 'Starting to obtain Syn/Nonsyn variant proportion per gene...'
    logging.info(t)
    # Define variants using environ
    environ['bed_variants'] = bed_variants
    environ['vcf_file'] = vcf_file
    environ['bed_genes'] = bed_genes
    environ['bed_intersect'] = bed_intersect
    # Generate bed_variants file
    l = 'bcftools query -f \'%CHROM\t%POS\t%END\t%INFO/ANN\n\' $vcf_file'
    l += ' | awk \'BEGIN{OFS=\"\t\"} {print $1, $2-1, $2, $4}\''
    l += ' > $bed_variants'
    run_silent(l, '')
    logging.info(f'{bed_variants} created.')
    l = f'bedtools intersect -a {bed_variants}'
    l += f' -b {bed_genes} -wa -wb > {bed_intersect}'
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
    # Initialize defaultdict (dictionary) for counts
    gene_counts = defaultdict(
        lambda: {"synonymous": 0, "nonsynonymous": 0}
        )
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
    return gene_counts

def extract_regions(bed_file:str) -> list[list[str,int,int,str]]:
    """
    Extracts regions of the genome from bed_file and formats them.
    Returns a list of regions with format chr, start, end, name.
    """
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
    return regions

def fragment_lengths(output_dir:str, sra_id:str, bam_sorted:str) -> tuple[float]:
    logging.info('Starting to obtain fragment lengths...')
    l = 'samtools view -f '
    l += f'0x2 {bam_sorted} |'
    l += ' awk \'{if ($9>0 && $9<1000) print $9}\''
    l += f' > {output_dir}/fl_{sra_id}.txt'
    #run_silent(l, log_file)
    os.system(l)
    # Open the created file to obtain values
    fragment_lengths = open(f"{output_dir}/fl_{sra_id}.txt", "r").read().splitlines()
    # Convert to integers
    fragment_lengths = list(map(int, fragment_lengths))
    # Get mean, median and standard deviation of fragment lengths (fl)
    fl_mean = statistics.mean(fragment_lengths)
    fl_median = statistics.median(fragment_lengths)
    fl_stdv = statistics.stdev(fragment_lengths)
    logging.info('Fragment mean, median and standard deviation:')
    logging.info(f'{fl_mean} / {fl_median} / {fl_stdv}')
    return fl_mean, fl_median, fl_stdv

def parse_ann_field(ann_field:str) -> list[str]:
    """
    Parses the ANN field and returns a list of effects.
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

def variants_per_bin(vcf_file:str, genome_sizes:str, bin_size:int) -> dict:
    # Initialize the output dictionary
    output_dict = {}
    # Generate genome bins
    genome = BedTool().window_maker(g=genome_sizes, w=bin_size)
    variants = vcf_to_bed(vcf_file)
    # Sort (just in case)
    genome = genome.sort()
    variants = variants.sort()
    # Intersect genome and variants to obtain variants per bin
    counts = genome.intersect(variants, c=True)
    if len(list(counts)) == 0:
        w = 'counts variable is empty. Variants per bin not recorded.'
        logging.warning(w)
        # Show the first 5 lines of genome bins
        print("GENOME BINS (BED format, first 5):")
        for i, feature in enumerate(genome):
            if i >= 5: break
            print(str(feature).strip())
        # Show the first 5 lines of the VCF input (converted to BED)
        print("\nVCF VARIANTS (BEDTools-parsed, first 5):")
        for i, feature in enumerate(variants):
            if i >= 5: break
            print(str(feature).strip())
        # Show the first 5 lines of the intersect result
        print("\nINTERSECT RESULT (first 5):")
        for i, feature in enumerate(counts):
            if i >= 5: break
            print(str(feature).strip())
    else:
        # Put counts in dataframe
        l_names = ['chrom', 'start', 'end', 'variant_count']
        df_vg_bins = counts.to_dataframe(names=l_names)
        if len(list(df_vg_bins)) != 0:
            df_vg_bins['bin_region'] = df_vg_bins.apply(
                lambda row: f'{row["chrom"]}:{row["start"]}-{row["end"]}',
                axis=1
                )
            total_varcount = sum(df_vg_bins['variant_count'])
            for row in df_vg_bins.itertuples():
                key = f'variants_in_{row.bin_region}'
                output_dict[key] = row.variant_count
                key = f'variants_in_{row.bin_region}_normalized'
                output_dict[key] = float(row.variant_count)/total_varcount
        else:
            w = 'df_vg_bins variable is empty.'
            w += ' Variants per bin not recorded.'
            logging.warning(w)
            # Show the first 5 lines of genome bins
            print("GENOME BINS (BED format, first 5):")
            for i, feature in enumerate(genome):
                if i >= 5: break
                print(str(feature).strip())
            # Show the first 5 lines of the VCF input (converted to BED)
            print("\nVCF VARIANTS (BEDTools-parsed, first 5):")
            for i, feature in enumerate(variants):
                if i >= 5: break
                print(str(feature).strip())
            # Show the first 5 lines of the intersect result
            print("\nINTERSECT RESULT (first 5):")
            for i, feature in enumerate(counts):
                if i >= 5: break
                print(str(feature).strip())
    return output_dict

def vcf_to_bed(vcf_file):
    def converter():
        for feature in BedTool(vcf_file):
            chrom = feature.chrom
            start = int(feature.start)
            ref = feature[3]
            end = start + max(len(ref), 1)
            yield f"{chrom}\t{start}\t{end}"
    return BedTool(converter()).saveas()
