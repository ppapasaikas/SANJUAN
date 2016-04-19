#!/bin/bash

mergedJunctions=$1
transcriptsBed=$2
output=$3

bedtools intersect -s -f 1 -wa -wb -a $mergedJunctions -b $transcriptsBed > $output

