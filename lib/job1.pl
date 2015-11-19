#!/usr/bin/perl

#$d6='\'\$6';
#$d1='\$1';
$d6='\'$6';
$d1='$1';
open(OUT,">".$ARGV[1]);
print OUT "samtools view -h $ARGV[0] | awk $d6 !~ /N/ || $d1 ~ /^@/\' | samtools view -bS - \n";
print OUT "sleep 300";
close(OUT);

# wait 5 min so that output on cluster is written completely for sure 
sleep 300;
