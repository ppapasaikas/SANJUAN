#!/usr/bin/perl

#$d6='\'\$6';
#$d1='\$1';
$d6='\'$6';
$d1='$1';
open(OUT,">".$ARGV[1]);
# the tr command should look in the end like this: perl -pe "'tr/\-\+/\+\-/ if $_=~/\/1\s/'"
print OUT "samtools view -h $ARGV[0] | awk $d6 !~ /N/ || $d1 ~ /^@/\' | samtools view -bS - | bamToBed -i - | perl -pe \'tr/\-\+/\+\-/ if \$_=~/\\/1\\s/\' - > $ARGV[2]\n";
print OUT "sleep 60";
close(OUT);

sleep 60;

chmod 0755, $ARGV[1];

# wait 5 min so that output on cluster is written completely for sure 
sleep 60;
