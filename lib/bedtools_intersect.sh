#!/bin/bash

rnaseq=$1  # S stranded U unstranded
out_intr_segm_sorted=$2
no_junct_reads=$3
output=$4
####-F: Fraction of read that needs to overlap.
####-f: Fraction of intron that needs to be covered by read.
####-e: Either F or f need be satisfied. Not necessarily both. Note: As of bedtools 2.25 this does not seem to work. i.e -e does not OR
if [ "$rnaseq" = "U" ]; then bedtools intersect -c -F 0.1 -e -split -sorted -a $out_intr_segm_sorted -b $no_junct_reads > $output; fi
if [ "$rnaseq" = "S" ]; then bedtools intersect -c -F 0.1 -e -split -s -sorted -a $out_intr_segm_sorted -b $no_junct_reads > $output; fi
