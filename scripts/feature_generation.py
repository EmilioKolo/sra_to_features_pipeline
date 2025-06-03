#!/usr/bin/env python3


from pybedtools import BedTool
import logging
import os
import statistics


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

def variants_per_bin(vcf_file:str, genome_sizes:str, bin_size:int) -> dict:
    # Initialize the output dictionary
    output_dict = {}
    # Generate genome bins
    genome = BedTool().window_maker(g=genome_sizes, w=bin_size)
    variants = BedTool(vcf_file)
    # Sort (just in case)
    genome = genome.sort()
    variants = variants.sort()
    # Intersect genome and variants to obtain variants per bin
    counts = genome.intersect(variants, c=True)
    if len(list(counts)) == 0:
        logging.warning('counts variable is empty. Variants per bin not recorded.')
    else:
        # Put counts in dataframe
        l_names = ['chrom', 'start', 'end', 'variant_count']
        df_vg_bins = counts.to_dataframe(names=l_names)

        if len(list(df_vg_bins)) == 0:
            df_vg_bins['bin_region'] = df_vg_bins.apply(lambda row: row['chrom'] + ':' + str(row['start']) + '-' + str(row['end']), axis=1)

            total_varcount = sum(df_vg_bins['variant_count'])

            for row in df_vg_bins.itertuples():
                key = f'variants_in_{row.bin_region}'
                output_dict[key] = row.variant_count
                key = f'variants_in_{row.bin_region}_normalized'
                output_dict[key] = float(row.variant_count)/total_varcount
        else:
            logging.warning('df_vg_bins variable is empty. Variants per bin not recorded.')
    return output_dict
