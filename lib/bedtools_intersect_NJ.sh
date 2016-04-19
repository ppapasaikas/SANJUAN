#!/bin/bash

mergedJunctions_pm10=$1
mergedJunctions_pm1000=$2
output=$3

bedtools intersect -s -wa -wb -a $mergedJunctions_pm10 -b $mergedJunctions_pm1000 > $output

