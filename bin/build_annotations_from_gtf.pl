#!/usr/bin/perl
use Cwd;
####### Builds SANJUAN required annotation files from gtf.
####### This only needs to be done in the case of adding a non-provided genome/genome_version
####### Usage: perl build_annotations_from_gtf.pl transcriptome_gtf prefix
### e.g: perl build_annotations_from_gtf.pl hg38.gtf hg38 

if (@ARGV<3 || $ARGV[0]=~ /help/) {
	print "Creates SANJUAN annotation files for new genomes\n";
	print "perl build_annotations_from_gtf.pl transcriptome_gtf prefix sanjuan_dir\n";
	print "transcriptome_gtf      : GTF files with gene annotation of genome to be added\n";
	print "prefix                 : short name of new genome; This name is used as ID to select this genome when you run SANJUAN later\n";
	print "sanjuan_dir            : full path to the installation directory of SANJUAN; This directory should have as sub-directory: bin, db, lib, perllib, mapping_indexes\n";
	exit 0; 
}

my $sanjuan_dir=$ARGV[2];
unless(-d "$ARGV[2]/db"){die "Cannot find sub-directory db in specified SANJUAN directory $sanjuan_dir\n";}

$OUTBED="$sanjuan_dir/db/SANJUAN_annotation_files/" . $ARGV[1] . "_Transcripts.bed";
$OUTT2J="$sanjuan_dir/db/SANJUAN_annotation_files/" . $ARGV[1] . "_Transcript_Junctions.txt";
$OUTT2I="$sanjuan_dir/db/SANJUAN_annotation_files/" . $ARGV[1] . "_TxID2Name.txt";
open (OUTBED, ">$OUTBED") || die;
open (OUTT2J, ">$OUTT2J") || die;
open (OUTT2I, ">$OUTT2I") || die;


####Build transcript bed file
print "\nBuilding transcript bed file...\n";

open (IN, $ARGV[0]) ||die;
	while (<IN>){
	$line=$_;
	chomp $line;
	@mat=split /\t/,$line;
	next unless $mat[2]=~/exon/i;
	$TxID=$1 if $line=~/transcript_id \"([^\"]+)\"/;
	$Gene{$TxID}=$1 if $line=~/gene_name \"([^\"]+)\"/;
	push @{$STARTS{$TxID}},$mat[3];
	push @{$ENDS{$TxID}},$mat[4];
	$CHR{$TxID}=$mat[0];
	$STRAND{$TxID}=$mat[6];
	}
close IN;


foreach $TxID (keys %STARTS) {
	$TxStart{$TxID}=min(@{$STARTS{$TxID}});
	$TxEnd{$TxID}=max(@{$ENDS{$TxID}});
	print OUTBED "$CHR{$TxID}\t$TxStart{$TxID}\t$TxEnd{$TxID}\t$TxID\t0\t$STRAND{$TxID}\n";	#Generate Bedfile
	}

print `sort -k1,1 -k2n,2 -k3n,3 -o $OUTBED $OUTBED`;
print "\nDone\n";


######Build transcript Junction file
print "\nBuilding transcript junction file...\n";
foreach $TxID (keys %STARTS) {
	@S = sort { $a <=> $b } @{$STARTS{$TxID}};
	@E = sort { $a <=> $b } @{$ENDS{$TxID}};
	for $i (0..$#S-1){
	  $end=$S[$i+1]-1;	#To match the used junction convention
	  print OUTT2J $TxID ."\t" . $CHR{$TxID} . '_' . $E[$i]  . '_' . $end  . '_' . $STRAND{$TxID} . "\n";	#Generate Tx2Jn table
	  }
	}
print `sort -k1,1 -o $OUTT2J $OUTT2J`;
print "\nDone\n";
close IN;


######Build transcript 2  GeneName File
print "\nBuilding transcript to gene name file...\n";
%trID2gN=();
$trID2gID=();
open (IN, $ARGV[0]) ||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$GeneName="";
$GeneID="";
next unless $mat[8]=~/transcript_id \"([^\"]+)\";/;
$TxID=$1;
next unless $mat[8]=~/gene_id \"([^\"]+)\";/;
$GeneID=$1;
$GeneName=$1 if $mat[8]=~/gene_name \"([^\"]+)\";/;


if(!defined($trID2gN{$TxID})){$trID2gN{$TxID}=$GeneName;}
if(!defined($trID2gID{$TxID})){$trID2gID{$TxID}=$GeneID;}
}

foreach $trID (keys %trID2gN){
	if($trID2gN{$trID}){print OUTT2I "$trID\t$trID2gN{$trID}\n";}
	else{
		# if we don't have a gene name, we take the gene id
		if($trID2gID{$trID}){print OUTT2I "$trID\t$trID2gID{$trID}\n";}
	}
}


print `sort -k1,1 -o $OUTT2I $OUTT2I`;
print "\nDone\n";


#######Retrieving genome file from UCSC
print "\nRetrieving genome file (chromosome sizes) from UCSC...\n";
$OUTGENOME="$sanjuan_dir/db/genomes/" . $ARGV[1] . ".genome";
print `$sanjuan_dir/lib/fetchChromSizes $ARGV[1]  > $OUTGENOME`;




print "\nDone\n";





sub min{
@vals=@_;
$min=1e20;
foreach (@vals){
$min=$_ if $_<$min;
}
return $min;
}


sub max{
@vals=@_;
$max=-1e20;
foreach (@vals){
$max=$_ if $_>$max; 
}
return $max;
}

