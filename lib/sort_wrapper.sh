#!/bin/bash

# $1 input_bed
# $2 bam
# $3 output_bed

bedtools sort -faidx <(samtools view -H $2 | grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}') -i $1 > $3