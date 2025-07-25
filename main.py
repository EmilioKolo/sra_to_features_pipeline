#!/usr/bin/env python3


"""
Main function to run the pipeline from SRA IDs or fastq files to 
a set of defined features.
"""

from scripts.arg_handler import get_variables
from scripts.pipeline import *
from scripts.quality_control import *
from scripts.log_scripts import *


def main():
    """
    Function to run when __name__=="__main__"
    """

    # Get variables using the get_variables function
    dict_var = get_variables()

    ### Display
    log_print(
        f"### Starting pipeline with {dict_var['sra_id']}...",
        level='info',
        log_file=dict_var['log_print']
    )
    ###

    # Check if sra was given as input
    if dict_var['sra_given']:
        # Download fastq files
        dict_var = download_sra(dict_var)
    
    # Align the fastq files to the reference sequence
    dict_var = align_to_reference(dict_var)

    # Perform variant calling and analysis
    dict_var = variant_call_and_analysis(dict_var)

    ### Display
    log_print(
        f'### Starting checks of intermediate files...',
        level='info',
        log_file=dict_var['log_print']
    )
    ###

    # Perform checks on the intermediate files
    perform_checks(dict_var)

    ### Display
    log_print(
        f'### Starting feature generation...',
        level='info',
        log_file=dict_var['log_print']
    )
    ###

    # Generate features
    dict_features = feature_generation(dict_var)

    ### Display
    log_print(
        f'### Features generated. Saving...',
        level='info',
        log_file=dict_var['log_print']
    )
    ###

    # Save features
    save_features(dict_var, dict_features)

    return None


if __name__=='__main__':
    main()
