#!/bin/bash


# $1 input_bed
# $2 bam
# $3 output_bed

#sort -k1,1 -k2,2n $1 > $2

samtools view -H $2 | grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}' > names_bam1_tempfile_87askjhmsjx6542.txt
bedtools sort -faidx names_bam1_tempfile_87askjhmsjx6542.txt -i $1 > $3 
rm names_bam1_tempfile_87askjhmsjx6542.txt

#sleep 1m
