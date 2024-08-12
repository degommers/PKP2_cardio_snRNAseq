#!/bin/bash

# AUTHOR: Demi Gommers
# Copyright: April 2024, UMC Utrecht, Utrecht University

# Stop execution of script whenever a non-zero exit status is returned
set -e

# Tools used: cellranger bamtofastq v1.4.1
# Script is run on 8core CPU, 32GB RAM

# Help function
Help() {
        # Display Help
        echo "This script converts bam files into fastq files with cellranger bamtofastq"
        echo "Required for this script are: directory with samples (sample_ID/Tissue/file.bam)"
	echo "Output directory is the same as input directory"
	echo
	echo "IMPORTANT: make sure that cellranger is installed and up-to-date"
        echo
        echo "Syntax: Input requirements [i|c|o]"
        echo "options:"
        echo " -i INPUT_DIR            Directory with all samples, structured like: sample_ID/tissue/file.bam"
}

# Parse command line arguments
while getopts "i:h:" flag; do
        case "${flag}" in
                i) input_dir=${OPTARG} ;;
                h) Help; exit ;;
                \?) Help; echo "Invalid option: ~$OPTARG"; exit 1 ;;
        esac
done

# Execute cellranger bamtofastq
for sample in $input_dir/*; do
	echo "Processing sample: $sample"
	for tissue in $sample/*; do
		echo "Processing tissue: $tissue"

		# Bam file variable
		bam_file=$(find $tissue -type f -name "*.bam")

		echo "Converting $bam_file to fastq"
		# Run cellranger bamtofastq
		cellranger bamtofastq --nthreads 12 --traceback $bam_file $tissue/fastq/
		echo "$bam_file converted to FASTQ file and saved here: $tissue/"
	done
done
