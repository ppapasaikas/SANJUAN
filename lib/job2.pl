#!/usr/bin/perl

$d6='\'$6';
open(OUT,">".$ARGV[1]);
#print `echo "samtools view $ARGV[0] | awk $d6 ~ /N/' " > $ARGV[1]`;
print OUT "samtools view $ARGV[0] | awk $d6 ~ /N/' > $ARGV[2]\n";
print OUT "sleep 60";
close(OUT);

sleep 60;
chmod 0755, $ARGV[1];

# wait 10 min so that output on cluster is written completely for sure 
sleep 60;
