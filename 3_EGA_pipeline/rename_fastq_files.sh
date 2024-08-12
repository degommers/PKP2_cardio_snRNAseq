#!/bin/bash

# Author: Demi Gommers
# Copyright: UMC Utrecht, Utrecht University
# May 2024

# Renaming the converted fastq files with sample_ID and tag

for dir in samples/BO*/*/fastq_140524_2/*/*
 do
   echo $dir
   sample_id=$(echo $dir | cut -d "/" -f 5 | sed -e 's/premrna_0_1_//')
   echo $sample_id
   rename "s/bamtofastq/$sample_id/" "$dir"
 done
