# AUTHOR: Demi Gommers
# Copyright: April 2024, UMC Utrecht, Utrecht University

# Tools used: cellranger bamtofastq v1.4.1
# Script is run on 8core CPU, 32GB RAM


# Stop execution of script whenever a non-zero exit status is returned
set -e

#Help function
Help() {
# Display Help
        echo "This script converts bam files into fastq files with cellranger bamtofastq"
        echo "Required for this script are: directory with samples (sample_ID/Tissue/file.bam)"
        echo "Output directory is the same as input directory"
        echo
        echo "IMPORTANT: make sure that cellranger is installed and up-to-date"
        echo
        echo "Syntax: script.sh -i INPUT_DIR"
        echo "options:"
        echo " -i INPUT_LIST Directory with all samples, structured like: sample_ID/tissue/file.bam"
}

# Function to convert BAM files to FASTQ
convert_to_fastq() {
        local bam_file="$1"
        local output_dir="$2"
	echo "Converting ${bam_file} to fastq"
        # Run cellranger bamtofastq
        cellranger bamtofastq --nthreads 16 --traceback "$bam_file" "$output_dir"
        echo "$bam_file converted to FASTQ file and saved in: $output_dir"

}

# Export function for parallel execution
export -f convert_to_fastq

# Parse command line arguments
while getopts "i:h" flag; do
        case "${flag}" in
                i) input_list=${OPTARG} ;;
                h) Help; exit ;;
                ?) Help; echo "Invalid option: -$OPTARG"; exit 1 ;;
        esac
done

# Execute cellranger bamtofastq in parallel
cat $input_list | parallel -j+0 convert_to_fastq {} '{= s|/[^/]+$|/fastq_140524|=}'
