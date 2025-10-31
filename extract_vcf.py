#!/usr/bin/env python3

"""
Script: extract_vcf.py

Description:
Decompresses a VCF.GZ file and extracts key information 
(CHROM, POS, REF, ALT, QUAL, ALT_DP, REF_DP) into a CSV file.

Usage:
    python3 extract_vcf.py <input_file.vcf.gz> <output_folder/>

Arguments:
    <input_file.vcf.gz>: Path to the compressed VCF file.
    <output_folder/>: Directory where the resulting CSV will be saved. 
                      The directory will be created if it does not exist.

Example:
    python3 extract_vcf.py data/sample.vcf.gz processed_output/
"""

import gzip
import csv
import os
import sys


def process_vcf_to_csv(vcf_gz_path, output_dir):
    """
    Decompresses a VCF.GZ file, extracts variant information, and saves it to a CSV file.

    Args:
        vcf_gz_path (str): Path to the input VCF.GZ file.
        output_dir (str): Path to the directory where the CSV should be saved.
    
    Returns:
        tuple: (success_status, message/path)
    """
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Define input/output filenames
    base_filename = os.path.basename(vcf_gz_path).replace('.vcf.gz', '')
    csv_filename = f"{base_filename}.csv"
    csv_output_path = os.path.join(output_dir, csv_filename)

    # The CSV columns to be extracted
    header = [
        'chromosome', 
        'variant position', 
        'reference nucleotide (REF)', 
        'variant nucleotide (ALT)', 
        'quality value (QUAL)', 
        'depth for the variant nucleotide (ALT_DP)', 
        'depth for non-variant nucleotides (REF_DP)'
    ]
    
    try:
        # Open the compressed VCF file for reading
        with gzip.open(vcf_gz_path, 'rt') as vcf_file:
            # Open the CSV file for writing
            with open(csv_output_path, 'w', newline='') as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerow(header)

                for line in vcf_file:
                    # Skip header lines, which start with '#'
                    if line.startswith('##'):
                        continue
                    
                    # The column header line starts with a single '#'
                    if line.startswith('#CHROM'):
                        # This line indicates we are at the data section. Skip it.
                        continue
                        
                    # Process variant data lines
                    try:
                        fields = line.strip().split('\t')
                        
                        # Extract basic fields
                        chrom = fields[0]
                        pos = fields[1]
                        ref = fields[3]
                        alt = fields[4]
                        qual = fields[5]
                        info_str = fields[7]
                        
                        # Parse INFO field for DP (Total Depth) and DP4 (Allelic Depths)
                        info_items = info_str.split(';')
                        info_dict = {}
                        for item in info_items:
                            if '=' in item:
                                key, value = item.split('=', 1)
                                info_dict[key] = value
                        
                        # DP4 format: N_ref_forward, N_ref_reverse, N_alt_forward, N_alt_reverse
                        dp4_str = info_dict.get('DP4')
                        
                        if dp4_str:
                            dp4_values = [int(x) for x in dp4_str.split(',')]
                            # Non-variant (Reference) depth = ref-forward + ref-reverse
                            ref_dp = dp4_values[0] + dp4_values[1]
                            # Variant (Alternate) depth = alt-forward + alt-reverse
                            alt_dp = dp4_values[2] + dp4_values[3]
                        else:
                            # Handle cases where DP4 might be missing (should be rare/filtered)
                            ref_dp = 'N/A'
                            alt_dp = 'N/A'
                            
                        # Write the extracted data to the CSV
                        csv_writer.writerow([
                            chrom, 
                            pos, 
                            ref, 
                            alt, 
                            qual, 
                            alt_dp, 
                            ref_dp
                        ])

                    except Exception as e:
                        print(f"Error processing line: {line.strip()}. Error: {e}")
                        continue
                        
        return True, csv_output_path

    except FileNotFoundError:
        return False, f"Error: Input VCF file not found at {vcf_gz_path}"
    except Exception as e:
        return False, f"An unexpected error occurred: {e}"


def main():
    # Check for correct number of arguments
    if len(sys.argv) != 3:
        print("Usage: python extract_vcf.py <input_file.vcf.gz> <output_folder/>", 
              file=sys.stderr)
        sys.exit(1)

    # Assign arguments
    input_vcf_path = sys.argv[1]
    output_folder = sys.argv[2]

    # Run the core logic
    success, result = process_vcf_to_csv(input_vcf_path, output_folder)
    
    if success:
        print(f"✅ Success! Variants extracted to: {result}")
        sys.exit(0)
    else:
        print(f"❌ Failure. {result}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
