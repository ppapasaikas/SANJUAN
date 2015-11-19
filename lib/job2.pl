#!/usr/bin/perl

$d6='\'$6';
open(OUT,">".$ARGV[1]);
#print `echo "samtools view $ARGV[0] | awk $d6 ~ /N/' " > $ARGV[1]`;
print OUT "samtools view $ARGV[0] | awk $d6 ~ /N/'\n";
print OUT "sleep 300";
close(OUT);

# wait 10 min so that output on cluster is written completely for sure 
sleep 300;
