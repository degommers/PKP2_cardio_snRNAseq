#!/bin/bash

# AUTHOR: Demi Gommers
# Copyright: May 2024, UMC Utrecht, Utrecht University

# IMPORTANT: This script on this VM only works when the proxies are turned OFF
# Turning off proxies: check /home/demi@mydre.org/Documents/commands_often_used.txt
# When proxies are turned off > white listed domains are not available

# This command runs in parallel

# Dependencies
# pyega3 v5.1.0, parallel, 

# Stop execution of script whenever a non-zero exit status is returned
set -e

# Help function
Help() {
        # Display Help
        echo "This script fetches data from the EGA Archive utilizing pyega3 tool based on your own credentials input and specific EGAF ID numbers"
        echo "Required for this script are: credentials_file.json, input file with EGAF IDs, and output directory"
        echo
        echo "Syntax: Input requirements [i|c|o]"
        echo "options:"
        echo " -i INPUT_FILE            File with EGAF IDs"
        echo " -o OUTPUT_DIR            Output directory that exists"
        echo " -c CREDENTIALS_FILE      Credentials file (.json) that has access to EGA Archive datasets"
}

# Parse command line arguments
while getopts "c:o:i:h:" flag; do
        case "${flag}" in
                c) credentials_file=${OPTARG} ;;
                i) input_file=${OPTARG} ;;
                o) output_dir=${OPTARG} ;;
                h)
                        Help
                        exit ;;
                \?)
                        Help
                        echo "Invalid option: ~$OPTARG"
                        exit 1 ;;
        esac
done

# Define fetch file function using pyega3
fetch_file() {
    local credentials_file="$1"
    local file="$2"
    local output_dir="$3"

    # Remove trailing carriage return characters
    file=$(echo ${file} | tr -d '\r')

    # Print file that is fetched
    echo "Fetching file: $file"

    # Executing pyega3 using 8 connections
    pyega3 -c 8 -cf "$credentials_file" fetch $file --output-dir "$output_dir"
}

# Export the fetch_file function and variables so they can be used by GNU Parallel
export -f fetch_file
export credentials_file
export output_dir

# Use GNU Parallel to run fetch_file in parallel for each file in input_file
cat "$input_file" | parallel --jobs 32 --no-notice fetch_file "$credentials_file" {} "$output_dir"
