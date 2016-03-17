#!/bin/bash

rnaseq=$1  # S stranded U unstranded
out_intr_segm_sorted=$2
no_junct_reads=$3
output=$4

if [ "$rnaseq" = "U" ]; then bedtools intersect -c -sorted -a $out_intr_segm_sorted -b $no_junct_reads > $output; fi
if [ "$rnaseq" = "S" ]; then bedtools intersect -c -s -sorted -a $out_intr_segm_sorted -b $no_junct_reads > $output; fi