#!/bin/bash

# Author: Demi Gommers

# Dependencies:
# FastQC

find samples/ -type f -name "*Y*.fastq.gz" | parallel -j+0 --bar 'echo {}; echo "Running fastqc"; fastqc --memory 10000 --threads 32 {}'
