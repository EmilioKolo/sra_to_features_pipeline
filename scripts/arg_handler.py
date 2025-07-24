#!/usr/bin/env python3


"""
Functions associated with defining I/O variables.
"""

import argparse
import configparser
import os


def define_bname(
        path_1:str,
        path_2:str,
        remove_folder:bool=True
    ) -> str:
    """
    Defines the basename between two file paths by getting the part they
    have in common. By default, removes the folder if it is included.
    """
    # Define if folder is removed
    if remove_folder:
        str_1 = path_1.split('/')[-1]
        str_2 = path_2.split('/')[-1]
    else:
        str_1 = path_1
        str_2 = path_2
    # Initialize the string that is returned
    bname = ''
    # Go through both strings
    for i in range(min(len(str_1), len(str_2))):
        if str_1[i]==str_2[i]:
            bname += str_1[i]
        else:
            break
    return bname

def getargs() -> argparse.Namespace:
    """
    Generates a parser and defines arguments. 
    Returns a Namespace element with the defined arguments.
    """
    # Initialize argument parser
    parser = argparse.ArgumentParser(
        description="Bioinformatics pipeline for SRA/FASTQ data."
    )
    # Define mutually exclusive arguments for SRA ID and fastq
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--sra-id", help="SRA Run Accession ID.")
    input_group.add_argument(
        "--fastq",
        nargs='+',
        help="One or two FASTQ file paths (single-end or paired-end)."
    )
    # Add an output dir argument
    h = "Directory to store output files"
    h += " (default: /content/data/output)."
    parser.add_argument(
        "--output-dir",
        default="/content/data/output",
        help=h
    )
    # Add number of threads as an argument
    h = "Number of threads to use for parallel processing"
    h += " (default: 1)."
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help=h
    )
    # Return parsed args
    return parser.parse_args()

def get_variables() -> dict[str:str|int|bool|list]:
    """
    Function that handles the definition of all needed variables.
    """
    # Initialize the dictionary that is returned
    dict_variables = {}
    # Define args and extract relevant arguments
    args = getargs()
    # Initialize configuration file parser
    config = configparser.ConfigParser()
    # Read the config.ini file
    config.read('config.ini')

    # Extract relevant variables from config.ini and args
    dict_variables['BASE_DIR'] = config['Paths']['BASE_DIR']
    dict_variables['fasta_ref'] = config['Paths']['REFERENCE_FASTA']
    dict_variables['gff_ref'] = config['Paths']['REFERENCE_GFF']
    dict_variables['kraken_db'] = config['Paths']['KRAKEN_DB']
    dict_variables['snpeff_dir'] = config['Paths']['SNPEFF_DIR']
    dict_variables['bed_file'] = config['Paths']['BED_FILE']
    dict_variables['bed_genes'] = config['Paths']['BED_GENES']
    dict_variables['genome_sizes'] = config['Paths']['GENOME_SIZES']
    dict_variables['genome_name'] = config['Parameters']['GENOME_NAME']
    dict_variables['bin_size_gvs'] = config['Parameters']['BIN_SIZE_GVS']
    dict_variables['THREADS'] = int(args.threads)
    dict_variables['OUTPUT_DIR'] = str(args.output_dir)
    # Define some base values
    dict_variables['DATA_DIR'] = os.path.join(
        dict_variables['BASE_DIR'],
        'data'
    )
    dict_variables['log_print'] = os.path.join(
        dict_variables['BASE_DIR'],
        'log/print.log'
    )
    dict_variables['log_scripts'] = os.path.join(
        dict_variables['BASE_DIR'],
        'log/scripts.log'
    )
    # Define if SRA ID or fastq were given
    if args.sra_id is not None:
        # SRA ID was provided
        dict_variables['sra_id'] = str(args.sra_id)
        dict_variables['sra_given'] = True
        dict_variables['l_fastq'] = []
        # Define a fastq directory
        dict_variables['fastq_dir'] = os.path.join(
            dict_variables['DATA_DIR'],
            'fastq'
        )
    elif args.fastq is not None:
        # FASTQ file(s) were provided
        dict_variables['sra_given'] = False
        # Check if one or more fastq files were given
        num_fastq_files = len(args.fastq)
        if num_fastq_files == 1:
            print("Single-end FASTQ file provided:", args.fastq[0])
            dict_variables['l_fastq'] = [str(args.fastq[0])]
            dict_variables['sra_id'] = dict_variables['l_fastq'][0].split(
                '/'
            )[-1].split('.')[0]
            # Define the fastq directory
            dict_variables['fastq_dir'] = os.path.dirname(args.fastq[0])
        elif num_fastq_files == 2:
            print(
                "Paired-end FASTQ files provided:",
                args.fastq[0],
                "and",
                args.fastq[1]
            )
            dict_variables['l_fastq'] = [
                str(args.fastq[0]),
                str(args.fastq[1])
            ]
            dict_variables['sra_id'] = define_bname(
                dict_variables['l_fastq'][0],
                dict_variables['l_fastq'][1]
            )
            # Define the fastq directories
            dir_fastq_0 = os.path.dirname(args.fastq[0])
            dir_fastq_1 = os.path.dirname(args.fastq[1])
            # Check that both fastq directories are the same
            if dir_fastq_1 != dir_fastq_0:
                er = 'The provided fastq files are in different folders.\n'
                er += dir_fastq_0 + '\n'
                er += dir_fastq_1
                raise ValueError(er)
            else:
                dict_variables['fastq_dir'] = dir_fastq_0
        else:
            er = "Expected 1 or 2 FASTQ files, but received "
            er += str(num_fastq_files)
            raise ValueError(er)
    
    # Define a list of fastq files with their directory included
    l_fastq_full = []
    for i in dict_variables['l_fastq']: # Empty if SRA ID is given
        l_fastq_full.append(os.path.join(dict_variables['fastq_dir'], i))
    dict_variables['l_fastq_full'] = l_fastq_full

    # Define non-base folders and files
    dict_variables['tmp_sra'] = os.path.join(
        dict_variables['DATA_DIR'],
        'tmp_'+dict_variables['sra_id']
    )
    dict_variables['sam_bam_dir'] = os.path.join(
        dict_variables['DATA_DIR'],
        'sam_bam'
    )
    dict_variables['vcf_dir'] = os.path.join(
        dict_variables['DATA_DIR'],
        'vcf'
    )
    dict_variables['counts_file'] = os.path.join(
        dict_variables['tmp_sra'],
        'counts.csv'
    )
    dict_variables['bed_variants'] = os.path.join(
        dict_variables['tmp_sra'],
        'variants.bed'
    )
    dict_variables['bed_intersect'] = os.path.join(
        dict_variables['tmp_sra'],
        'intersect.bed'
    )
    return dict_variables
