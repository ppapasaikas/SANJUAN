#!/usr/bin/perl
use strict;
use warnings;
################ Merge Junction Files

my $Tophat_Junctions1 = $ARGV[0];
my $Tophat_Junctions2 = $ARGV[1];
my $out_file=$ARGV[2];

my %count=();
my @fs;
my $fh;

open ($fh,"<".$Tophat_Junctions1);
while (<$fh>){
	chomp;
	@fs=split("\t");
	$count{$fs[3]}=$fs[4];
}
close($fh);
	
open ($fh,"<".$Tophat_Junctions2);
while (<$fh>){
	chomp;
	@fs=split("\t");
	$count{$fs[3]}+=$fs[4];
}
close($fh);

open ($fh,">".$out_file);
foreach (keys %count){
	@fs=split("_");
	print $fh "$fs[0]\t$fs[1]\t$fs[2]\t$_\t$count{$_}\t$fs[3]\n";
	}
close($fh);


# wait 5 min so that output on cluster is written completely for sure 
#sleep 300;
	