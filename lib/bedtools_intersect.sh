str_spec=$1
out_intr_segm_sorted=$2
no_junct_reads=$3
output=$4

bedtools intersect -c $str_spec -sorted -a $out_intr_segm_sorted -b $no_junct_reads > $output
