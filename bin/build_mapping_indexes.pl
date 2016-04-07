use strict;
use Cwd;
####### Builds bowtie and tophat indexes. Usage: perl build_mapping_indexes.pl genome_fasta prefix transcriptome_gtf
### e.g: perl build_mapping_indexes.pl hg38.fa hg38 

my $pwd=cwd();
unless ($pwd=~/SANJUAN\/bin.*$/){
	die "\nPlease run this script from within the SANJUAN/bin directory";
	}


if ($#ARGV<2) {
	die "\n\n!!Missing arguments!!. Usage:
	perl build_mapping_indexes.pl genome_fasta prefix transcriptome_gtf\n
	e.g: 
	perl build_mapping_indexes.pl hg38.fa hg38\n\n\n"; 
	}

my $fasta=$ARGV[0];
my $prefix=$ARGV[1];
my $gtf=$ARGV[2];
my $tophat_prefix=$prefix . '_transcripts';

print "\nCreating bowtie2 indexes in SANJUAN/db/mapping_indexes/$prefix/ with prefix $prefix ...\n";
print `mkdir -p db/mapping_indexes/$prefix/$prefix`;
print `bowtie2-build $fasta ../db/mapping_indexes/$prefix/$prefix`;
print "Done ...\n";
print "\nCreating tophat2 transcriptome indexes in SANJUAN/db/mapping_indexes/$prefix/ with prefix $tophat_prefix  ...\n";
print `tophat2 -G $gtf --transcriptome-index=../db/mapping_indexes/$prefix/$tophat_prefix ../db/mapping_indexes/$prefix/$prefix`;
print "Done ...\n";

