#!/bin/bash

# AUTHOR: Demi Gommers
# Copyright: May 2024, UMC Utrecht, Utrecht University


# This script renames the fastq files after conversion from bam to fastq, filename is now bamtofastq_S1_L001_I1_001.fa.gz
# Files are renamed by replacing 'bamtofastq' in the filename with the sample_id, tissue_type and tag

# Stop execution of script whenever a non-zero exit status is returned
set -e

# Help function
Help() {
        # Display Help
        echo "This script creates directories per sample_ID and moves the overlapping sample_ID to the right directory"
        echo "Required for this script: output directory and input directory"
	echo "Input directory should be the output of the pyega fetch command"
        echo
        echo "Syntax: Input requirements [o|i]"
        echo "options:"
        echo " -o OUTPUT_DIR            Output directory"
	echo " -i INPUT_DIR		Input directory that exists"
}

# Parse command line arguments
while getopts "o:i:h:" flag; do
        case "${flag}" in
		i) input_dir=${OPTARG} ;;
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


# Creating a directory structure based on the filenames
for dir in "$input_dir"/EGAF*/; do
	echo "Processing directory: $dir"

	# Looping through each file in the directory
	for file in $dir/*; do
		filename=$(basename $file)
		echo "Processing file: $filename"

		# Extract sample ID from file name as directory name
		directory=$(echo "$filename" | grep -oE 'BO_H[0-9]{2}')
		echo "Sample ID: $directory"

		# Extract tissue type from file name as subdirectory name
		subdirectory=$(echo "$filename" | grep -oE '(LV(0|1)|RV(0|1)|S0(0|N))')
		echo "Tissue type: $subdirectory"

		# Create directory with sample_id as name
		mkdir -p "$output_dir/$directory/$subdirectory"

		# Move file to its corresponding subdirectory
		mv "$dir/$filename" "$output_dir/$directory/$subdirectory/"
		echo "Moved file: $filename to $output_dir/$directory/$subdirectory"
	done
done
