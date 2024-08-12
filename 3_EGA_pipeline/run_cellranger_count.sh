#!/bin/bash

# Author: Demi Gommers


# Define the path to the parent directory containing all sample directories
parent_dir="/mnt/data/demi/snRNAseq_PKP2_control_samples/samples"

# Output CSV file
csv_file="cellranger_config_PKP2_Control_snRNAseq.csv"

# Create header for CSV file
# echo "#[libraries]" > "$csv_file"
echo "sample,fastqs,library_type" > "$csv_file"

# Iterate over parent directories
for sample_dir in "$parent_dir"/*; do
    if [ -d "$sample_dir" ]; then
        # Iterate over sample subdirectories
        for sub_dir in "$sample_dir"/*; do
	    for sub_sub_dir in "$sub_dir"/*; do
		for sub_sub_sub_dir in $sub_sub_dir/*; do
			if [ -d "$sub_sub_sub_dir" ]; then
        		  # Construct the path to the FASTQ directory for the sample
                	  fastq_dir="$sub_sub_sub_dir/"
			  # Extract sample name
			  echo "Creating config file"
			  echo "--- Extracting sample info ---"
			  sample=$(echo "$fastq_dir" | grep -oE 'BO_H[0-9]{2}' | head -n 1)
			  sample_origin=$(echo "$fastq_dir" | grep -oE '(LV(0|1)|RV(0|1)|S0(0|N))'| head -n 1)
			  sample_tag=$(echo "$fastq_dir" | cut -d "_" -f 11 | tr -d "/")
			  sample_name="${sample}_${sample_origin}_${sample_tag}"

			  echo "${sample_name} directory included"
                	  # Append sample information to the CSV file
                	  echo "$sample_name,$fastq_dir,Gene Expression" >> "$csv_file"
			  echo "Created CSV config file"
            		fi
		done
	    done
        done
    fi
done

# Path to transcriptome file
transcriptome_ref="/home/demi@mydre.org/Downloads/reference_genome/hg38/release_111/refdata_cellranger_8_GRCh38_111_premrna"

# Run Cell Ranger Count with the generated CSV file
echo "Running cellranger count with created ${csv_file}"
cellranger count --id="snRNAseq_count_cellranger_8_PKP2_Control" \
  --transcriptome="$transcriptome_ref" \
  --libraries="$csv_file" \
  --create-bam=false \
  --localcores=32 \
  --localmem=120 \
  --disable-ui 
