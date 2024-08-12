#!/bin/bash

# AUTHOR: Demi Gommers
# Copyright: April 2024, UMC Utrecht, Utrecht University

# IMPORTANT: This script on this VM only works when the proxies are turned OFF
# Turning off proxies: check /home/demi@mydre.org/Documents/commands_often_used.txt
# When proxies are turned off > white listed domains are not available

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
	echo " -i INPUT_FILE		File with EGAF IDs"
	echo " -o OUTPUT_DIR		Output directory that exists"
	echo " -c CREDENTIALS_FILE	Credentials file (.json) that has access to EGA Archive datasets"
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


# Define fetch files function using pyega3
fetch_files() {
	# Setting input variables for the function fetch_files
	local credentials_file="$1"
	local input_file="$2"
	local output_dir="$3"

	# While loop to execute pyega3 for every file in input_file
	while IFS= read -r file; do
		# Remove trailing carriage return characters
		file=$(echo ${file} | tr -d '\r')

		# Print file that is fetched
		echo "Fetching file: $file"
		
		# Executing pyega3 using 5 connections
		pyega3 -c 8 -cf "$credentials_file" fetch $file --output-dir "$output_dir"
	done < "$input_file"

}

# Call function with said arguments
fetch_files "$credentials_file" "$input_file" "$output_dir"
