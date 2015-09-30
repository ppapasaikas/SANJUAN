$genome="hg"; 	#Specify species genome: hg -> human, mm-> mouse, dr-> zebrafish

#Specify directory for Input:
$indir='/users/sdelaluna/sequencing_data/Chiara_DiVona/2014-08-27/';	# !!! Do not forget the last slash (/) !!!

#Specify base directory for Output:
$basedir='/users/jvalcarcel/ppapasaikas/CHIARA/';	# !!! Do not forget the last slash (/) !!!

$map_no_skip=0; #Go to mapping directly (i.e map using untrimmed fastq files)
$skip{1}=0;	#Skip individual step(s). Set this to 1 to skip trimming (i.e if trimming has already completed)
$skip{2}=0;	#Skip individual step(s). Set this to 1 to skip mapping (e.g to only perform trimming)






##########################################################################################################
##########################################################################################################
##########################################################################################################

@files = <$indir*>;
foreach (@files){
$f=$_;
#print "$f\n";
next unless (-f $f);
#push @READ1,$f if  $f=~/read1.+.gz/;
#push @READ2,$f if  $f=~/read2.+.gz/;
push @READ1,$f if  $f=~/1.fastq.gz$/;
push @READ2,$f if  $f=~/2.fastq.gz$/;
#push @READ1,$f if  $f=~/_R1*fq.gz$/;	#Manu's Files
#push @READ2,$f if  $f=~/_R2*fq.gz$/;	#Manu's Files
}

@READ1=sort(@READ1);
@READ2=sort(@READ2);

$nc=$#READ1+1;
print "\nFound $nc pair-end datasets:\n";
for $cf (0..$#READ1){
push @COND, $1 if $READ1[$cf] =~/$indir(.+)_\d{4,}/;
push @PREFIX, $1 if $READ1[$cf] =~/$indir(.+)_read/;

#push @COND, $1 if $READ1[$cf] =~/$indir(.+)_R[12]/;	#Manu's Files
#push @PREFIX, $1 if $READ1[$cf] =~/$indir(.+)_R[12]/;	#Manu's Files

##push @COND, $1 if $READ1[$cf] =~/$indir(.+)_[12].fastq/;
##push @PREFIX, $1 if $READ1[$cf] =~/$indir(.+)_[12].fastq/;
print "$1\n";
}


########################################	T R I M M I N G		##########################################
#Specify adapter sequence to be trimmed. To find adapter: minion search-adapter -i FASTQFILE.gz
$ad="AGATCGGAAGAGC";	#(Default adapter sequence for CRG facility, RNAseq )
$bc="";
$ad_bc=$ad . $bc;

$odir=$basedir . 'TRIM';

$skip{1}=1 if $map_no_skip==1;

unless ($skip{1}){
print "\nTrimming...\n";
unless (-d $odir){
print `mkdir $odir`;
}


for $cf (0..$#READ1){
$cond=$COND[$cf];
$jname=$cond;
$jname=~ s/^\d+//;
print `qsub -q long-sl65 -V -cwd  -N $jname._TRIM -e ~/temp -pe ompi 1 -l virtual_free=40G -l h_rt=48:00:00 -b y trim_galore --stringency 3 -q 0 -a $ad_bc -a2 $ad_bc --length 19 --paired $READ1[$cf] $READ2[$cf] -o $odir > flush.tmp`;
}
print `rm -rf *_TRIM*`;
}


########################################	M A P P I N G		##########################################
##TO GET BOWTIE INDEXES:
#print `wget -P /Volumes/HD2/RESEARCH/DATA/BOWTIE2_INDEXES/mm10/ ftp://ftp.cbcb.umd.edu/pub/data/bowtie2_indexes/incl/mm10.zip`;	#hg19.zip/mm10.zip
#OR from ilumina's iGENOMES
#print `wget -P /Volumes/HD2/RESEARCH/DATA/BOWTIE2_INDEXES/mm10/ ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz`;
#then: tar -zxvf ....gz;

$fqdir=$basedir . 'TRIM/';
$fqdir=$indir if $map_no_skip==1;


unless ($skip{2}){
print "\nMapping...\n";
for (0..$#READ1){
$cond=$COND[$_];
$jname=$cond;
$jname=~ s/^\d+//;
$fq1=$fqdir . $PREFIX[$_] . '_read1_val_1.fq.gz';
$fq2=$fqdir . $PREFIX[$_] . '_read2_val_2.fq.gz';
#$fq1=$fqdir . $PREFIX[$_] . '_read1.fastq.gz' if $map_no_skip==1;
#$fq2=$fqdir . $PREFIX[$_] . '_read2.fastq.gz' if $map_no_skip==1;
$fq1=$READ1[$_] if $map_no_skip==1;
$fq2=$READ2[$_] if $map_no_skip==1;


$top_out=$basedir . 'TOPHAT_' . $cond;
$hold=$jname .'._TRIM' ;
print "$top_out\n";


#Default Human:
print `qsub -q long-sl65 -V -cwd -N $jname._TOP -hold_jid $hold -e ~/temp -pe ompi 12 -l virtual_free=64G -l h_rt=48:00:00 -b y tophat2 --no-mixed --library-type fr-firststrand -o $top_out --transcriptome-index=/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/ -G /users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/cuffcmp.combined.gtf -p 12 -r 25 --mate-std-dev 85 -i 50 -I 800000 -x 1 /users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/hg19/hg19 $fq1 $fq2` if $genome eq "hg";

#Default Mouse:
print `qsub -q long-sl65 -V -cwd -N $jname._TOP -hold_jid $hold -e ~/temp -pe ompi 12 -l virtual_free=64G -l h_rt=48:00:00 -b y tophat2 --no-mixed --library-type fr-firststrand -o $top_out --transcriptome-index=/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/ENSEMBL_mm10_GRVm30/ -G /users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/ENSEMBL_mm10_GRVm30/cuffcmp.combined.gtf -p 12 -r 25 --mate-std-dev 85 -i 50 -I 800000 -x 1 /users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/mm10/mm10 $fq1 $fq2` if $genome eq "mm";

#Default Zebrafish:
print `qsub -q long-sl65 -V -cwd -N $jname._TOP -hold_jid $hold -e ~/temp -pe ompi 12 -l virtual_free=64G -l h_rt=48:00:00 -b y tophat2 --no-mixed --library-type fr-firststrand -o $top_out --transcriptome-index=/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/danRer10_ENSEMBL/ -G /users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/danRer10_ENSEMBL/cuffcmp.combined.gtf -p 12 -r 25 --mate-std-dev 85 -i 50 -I 800000 -x 1 /users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/dr10/dr10 $fq1 $fq2` if $genome eq "dr";

#Only for CNAG sequencing
##print `qsub -q long-sl65 -V -cwd -N $cond._TOP -hold_jid $hold -e ~/temp -pe ompi 12 -l virtual_free=64G -l h_rt=48:00:00 -b y tophat2 --no-mixed --library-type fr-unstranded -o $top_out --transcriptome-index=/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/ -G /users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/cuffcmp.combined.gtf -p 12 -r 20 --mate-std-dev 45 -i 50 -I 800000 -x 1 /users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/hg19/hg19 $fq1 $fq2`;

}
print `rm -rf *_TOP*`;
}











