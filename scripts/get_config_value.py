#!/usr/bin/env python3


"""
Function to be able to use config.ini variables for bash scripts.
"""

import configparser
import sys


# Require a strict number of input arguments
if len(sys.argv) != 4:
    print(
        "Usage: get_config_value.py <config_file> <section> <key>",
        file=sys.stderr
    )
    sys.exit(1)

# Define config file, section and key
config_file = sys.argv[1]
section = sys.argv[2]
key = sys.argv[3]

# Start the configparser parser
config = configparser.ConfigParser()
try:
    # Open config file
    config.read(config_file)
    # Check if the file has the corresponding section and key
    if config.has_section(section) and config.has_option(section, key):
        print(config.get(section, key))
    else:
        print(
            'The following key and section could not be found:',
            f'key: {key}\nsection: {section}',
            file=sys.stderr
        )
        sys.exit(1)
except Exception as e:
    print(f"Error reading config file: {e}", file=sys.stderr)
    sys.exit(1)
