$dir_out=$ARGV[0];
$dir_sj=$ARGV[1];


$output=`perl $dir_sj/consolidate_results.pl $dir_out/Annotated_Diff_Junctions.txt`;
open($out,">$dir_out/Annotated_Diff_Junctions_clean.txt") or die;
print $out $output;
close($out);

$output=`perl $dir_sj/consolidate_results.pl $dir_out/Annotated_Diff_Junctions_NC.txt`;
open($out,">$dir_out/Annotated_Diff_Junctions_NC_clean.txt") or die;
print $out $output;
close($out);

$output=`perl $dir_sj/SANJUAN_analyze_CEx.pl $dir_out/Annotated_Diff_Junctions_clean.txt $dir_out/Diff_Junctions_NC.txt`;
open($out,">$dir_out/CEs_clean.txt") or die;
print $out $output;
close($out);

$output=`perl $dir_sj/SANJUAN_analyze_CEx.pl $dir_out/Annotated_Diff_Junctions_NC_clean.txt $dir_out/Diff_Junctions_NC.txt`;
open($out,">$dir_out/CEs_NC_clean.txt") or die;
print $out $output;
close($out);                                                                                                             

# replace all entries in CEs_NC_clean.txt with entries in (more strict) CEs_clean.txt
%ids=();
open($fh,"$dir_out/CEs_NC_clean.txt") or die;
$header=<$fh>;
while(<$fh>){@fs=split("\t");$ids{join("_",@fs[1..4])}=$_;}
close($fh);

open($fh,"$dir_out/CEs_clean.txt") or die;
$header=<$fh>;
while(<$fh>){@fs=split("\t");$ids{join("_",@fs[1..4])}=$_;}  # add or overwrite if already existent 
close($fh);


# rewrite file
open($fh,">$dir_out/CEs_NC_clean.txt") or die;
foreach $key (keys %ids){print $fh $ids{$key};}
close($fh);