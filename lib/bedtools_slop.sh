input_bed=$1
genome_path=$2
b_parameter=$3
output_bed=$4

bedtools slop -i $input_bed -g $genome_path -b $b_parameter > $output_bed