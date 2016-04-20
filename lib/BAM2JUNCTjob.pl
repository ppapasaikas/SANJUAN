#!/usr/bin/perl
use strict;
use warnings;
my $sanjuan_dir=$ARGV[0];
my $OUTsh=$ARGV[1];
my $bamfile=$ARGV[2];
my $low_seq_req=$ARGV[3];
my $JunctFile=$ARGV[4];


my $d6='\'$6';
open(OUT,">".$OUTsh);
print OUT "samtools view $bamfile | awk $d6 ~ /N/'| perl $sanjuan_dir/get_juncts.pl $low_seq_req $JunctFile\n";
print OUT "sleep 60";
close(OUT);
sleep 10;
chmod 0755, $ARGV[1];
sleep 10;
