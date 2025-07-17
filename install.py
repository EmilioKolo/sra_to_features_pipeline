#!/usr/bin/env python3


"""
Functions to be used to set up the Docker environment.
"""

import configparser


def main():
    # Initialize configuration file parser
    config = configparser.ConfigParser()
    # Read the config.ini file
    config.read('config.ini')
    # Define base variables
    genome_sizes = config['Paths']['GENOME_SIZES']
    fasta_file = config['Paths']['REFERENCE_FASTA']
    gff_file = config['Paths']['REFERENCE_GFF']
    bed_genes = config['Paths']['BED_GENES']

    # Creation of genome.sizes file
    print(f'Creating {genome_sizes} file...')
    with open(fasta_file, 'r') as f:
        dict_fasta = {}
        name = None
        length = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save records when finding a new one
                if name is not None:
                    dict_fasta[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
        # Save last record
        if name is not None:
            dict_fasta[name] = length
    # Save dict to file
    with open(genome_sizes, 'w') as f:
        for key, value in dict_fasta.items():
            str_out = f'{key}\t{value}\n'
            f.write(str_out)

    print(f'{genome_sizes} file created.')

    # Create list of genes from gff file
    print(f'Extracting gene_data from gff file.')
    gene_data = []
    with open(gff_file) as f:
        # Go through gff file
        for line in f:
            # Skip header/comments
            if line.startswith("#"):
                continue
            # Turn line into list
            parts = line.strip().split("\t")
            # Skip lines that don't fit the expected format
            if len(parts) != 9:
                continue
            # Define feature type
            feature_type = parts[2]
            # Select genes
            if feature_type != "gene":
                continue
            # Define chromosome, start and end positions
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            # Select attributes
            attributes = parts[8]
            # Extract gene name
            gene_name = None
            for key in ["gene_name", "Name", "ID"]:
                prefix = key + "="
                if prefix in attributes:
                    # Get the first matching value
                    gene_name = attributes.split(prefix)[1].split(";")[0]
                    break
            # If there is a gene name, append it to gene_data
            if gene_name:
                gene_data.append([chrom, start, end, gene_name])

    print(f'Creating {bed_genes} file...')
    # Overwrite bed_genes as an empty file
    with open(bed_genes, 'w') as f:
        pass
    # Reopen in append mode
    with open(bed_genes, 'a+') as f:
        # Save genes from gene_data
        for gene in gene_data:
            chr = gene[0]
            start = gene[1]
            end = gene[2]
            name = gene[3]
            str_bed = f"{chr}\t{start}\t{end}\t{name}\n"
            f.write(str_bed)
    print(f'{bed_genes} file created and populated.')

if __name__=='__main__':
    main()
