#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path cwd);

### Verion
my $version="1.0 beta";

#### main paths parameters
# is the full path to the script sanjuan.pl including sanjuan.pl
# even, if sanjuan was called through a link 
my $abs_path=abs_path($0);
# remove /bin/sanjuan.pl to obtain root directory of SANJUAN installation
$abs_path =~ s/\/bin\/sanjuan\.pl//;
my $sanjuan_dir=$abs_path."/lib";   			# should contain all scripts used by SANJUAN
my $sanjuan_perllib=$abs_path."/perllib";		# should contain Perl modules Text and Statistics XXX check if we need both!
my $sanjuan_genomic_data_dir=$abs_path."/db";	# should contain folders genomes and annotation_files


# Are all important parts of SANJUAN in place?
# check if all sanjuan files can be accessed
my @SANJUAN_files=("annotate_Diff_Used_Introns.pl","annotate_Diff_Used_Junctions.pl","BAM2JUNCTjob.pl","bedtools_intersect_NJ.sh","bedtools_intersect_ST.sh","bedtools_intersect.sh","bedtools_slop.sh","calc_INTRON_retention.pl","calc_JUNCT_efficiency.pl","fetchChromSizes","get_juncts.pl","ln_s_wrapper.sh","Mapped_junctions2IntronSegments.pl","merge_junctions.pl","preProcess_and_Map.pl","SANJUAN_wrapper.pl","sort_wrapper.sh");
foreach (@SANJUAN_files){unless(-r $sanjuan_dir."/".$_){die "Cannot open file $_ which is essential for SANJUAN. You might go through the installation process again to solve this problem.\n";}}

# check if sub-directories genomes and annotation_files exist in db sub-directory
if(!-d $sanjuan_genomic_data_dir."/genomes" || !-d $sanjuan_genomic_data_dir."/SANJUAN_annotation_files"){die "Sub-directories genomes and/or annotation_files cannot be found in $sanjuan_genomic_data_dir. Try to specify parameter DBLOCATION / -db to specify their location\n";}

# check if sub-directories perllib is in place
if(!-d $sanjuan_perllib."/Text" || !-d $sanjuan_perllib."/Statistics"){die "Sud-directories / perl modules Text and/or Statistics cannot be found in $sanjuan_perllib. You might go through the installation process again to solve this problem.\n";}
#############################################################
#############################################################


sub print_help{
	print "\nSANJUAN $version -- *S*plicing *AN*alysis & *JU*nction *AN*notation\n================\n\n";
	print "Parameters can be defined by a parameter file or via command line arguments.\n\n";
	print "* parameters via parameter file, run\n\tsanjuan <parameter-file>\n";
	print "\tTo create a parameter file with explanations (file SANJUAN_parameters.txt will be created) run: sanjuan -exampleF\n";
	print "\n";
	print "* parameters via command line arguments, run\n";
	print "\tsanjuan ... ARGUMENTS ...\n\n";	#CinS
	print "\tARGUMENTS:\n";
	print "\t-nprocs number of parallel processes\n"; 
	print "\t-g     genome (ID); To see installed genomes run sanjuan -g.\n";
	print "\t-b     M -> start with mapping, S -> start with splicing analysis SKIPPING mapping\n";	#CinS
	print "\t-g1    short name for group 1\n";
	print "\t-f1    input files for group 1; Depending on value of argument -b, -f1 defines different types of input files.\n";
	print "\t\t\t if -b is M: pairs of FASTQ/FASTQ.GZ/FASTQ.BZ2 files describing RNAseq data,\n";
	print "\t\t\t\t comma separated list without white spaces following this order (file names don't matter): run1_read1.fastq,run1_read2.fastq,run2_read1.fastq,run2_read2.fastq,..\n";
	print "\t\t\t if -b set to S: exactly one BAM file containing all mapped reads for group 1\n";	#CinS
	print "\t-g2    like argument -g1 but for group 2\n";
	print "\t-f2    like argument -f1 but for group 2\n";
	print "\t-d     directoy of input data; If specified, files given via -f1 / -f2 are taken from this directory.\n";
	print "\t-o     output directory; If omitted current working directory is taken.\n";
	print "\t-p     encoding of base qualities in FASTQ files; values: phred33 (ASCII+33), phred64 (ASCII+64); Can be omitted if -b S.\n";
	print "\t\t\t For details have a look at section Encoding of Wikipedia article on FASTQ https://en.wikipedia.org/wiki/Fastq\n";	
	print "\t-l     RNAseq library type; values: 1S (single-end, stranded), 1U (single-end, unstranded), 2U (paired-end, unstranded), 2S (paired-end, stranded)\n";
   	print "\t\t\t Paired-end reads are expected in be given in two FASTQ files.\n"; #CinS  
   	print "\t\t\t Can be omitted if -b S.\n";
	print "\t-a     RNAseq adapter sequence; CRG standard is AGATCGGAAGAGC. Can be omitted if -b S.\n"; 
	print "\t\t\t To identify adapter you could try minion search-adapter -i FASTQFILE.gz\n";
	print "\t-tpm   If given it activates the Basic twopassMode STAR mapping option. More sensitive novel junction discovery at the cost of speed.\n"; #NinS
	print "\t-c     threshold on reported differentially spliced junctions; values VHC -> very high confidence (DPSI>20%, p-val<0.0001),\n";
	print "\t\t\t HC -> high confidence (DPSI>15%, p-val<0.001), MC -> medium confidence (DPSI>10%, p-val<0.01)\n";
	print "\t-i     High sensitivity intron retention analysis (IRM mode) will be done.\n";
	print "\t-s     Supporting junction evidence for IR identification (for IRM mode -i) will be required.\n";
	print "\t-lsr   Reads will be filtered-out less strictly (low sequence requirements).\n";
	print "\t\t\t If result files are empty, a reason could be too strong filtering. Then try to use -lsr.\n";
	print "\t-rmdup remove PCR duplicates\n";
	print "\t-t     program calls will be printed but not excecuted (test run).\n";
	print "\t-db    Full path to the directory db where SANJUAN will find pre-defined exon-exon junctions, genomes, and annotations.\n";
	print "\t\t\t Has to be set only if this directory is not under the SANJUAN installation directory.\n";
	print "\t-noqsub stand-alone run without sending jobs to CRG cluster\n\n";
	print "\tSTANDARD VALUES ALREADY SET:\n";
	print "\t-nprocs 1 -g hg19 -g1 COND -g2 CNTR -b M -o . -p phred33 -l 2S -a AGATCGGAAGAGC -c HC\n";
	print "\nFor printing a full example SANJUAN call: sanjuan -exampleC\n\n";
	print "Example call for human RNAseq data from CRG:\n";
	print "\tsanjuan -g1 ko -g2 cntr -f1 run1_1.fastq,run1_2.fastq -f2 run2_1.fastq,run2_2.fastq -d . -c HC -i -s\n\n";
	print "Contact:\n";
	print "Panagiotis Papasaikas started and developed SANJUAN: panagiotis.papasaikas\@crg.eu\n";
	print "Andre Gohr: andre.gohr\@crg.eu (Functionality Extensions and Support)\n\n";
}

if(@ARGV==0 || $ARGV[0] eq "--help" || $ARGV[0] eq "-help" || $ARGV[0] eq "help" || $ARGV[0] eq "?"){
	print_help;
	exit(0);
}

if(@ARGV==1 && $ARGV[0] eq "-exampleC"){
	print "\nsanjuan -o . -g hg19 -g1 COND -f1 <f1_read1,f1_read2,f2_read1,f2_read2..> -g2 CNTR -f2 <f1_read1,f1_read2,..> -p phred33 -l 2S -a AGATCGGAAGAGC -b M -c HC -i -s -lsr -noqsub\n\n";
	exit 0;
}

if(@ARGV==1 && $ARGV[0] eq "-exampleF"){
	if(-f "SANJUAN_parameters.txt"){print "File SANJUAN_parameters.txt already exists. Will not write example parameter file.\n";
		exit 0;
	}
	open(my $fh,">"."SANJUAN_parameters.txt") or die "Cannot open file SANJUAN_parameters.txt for writing";

	$\="\n";	
	print $fh "#######################  General notes on behavior of pipeline  #######################";
	print $fh "# 1. SANJUAN has two starting points from which you can start the analysis depending on the data you are working with."; #CinS 
	print $fh "# 1.a raw FASTQ files: specify the parameter RAWFASTQS in section Data to run mapping, and splicing analysis.";
	print $fh "# 1.b BAM files: specify the parameters BAM1 and BAM2 in section Data to SKIP trimming and mapping, AND ONLY RUN splicing analysis.";
	print $fh "# 2. SANJUAN will always start at the latest possible starting point. E.g., if you specify BAM1 and BAM2, SANJUAN will start directly with the splicing analysis regardless of whether you have specified RAWFASTQS or not.";	#CinS
	print $fh "#";
	print $fh "# 3. To avoid un-necessary computations, SANJUAN applies the following rules to decide which computations are to be skipped:";
	print $fh "# 4. If intermediate output files / results of preprocessing or splicing analysis are present, they will be used and not re-computed.";
	print $fh "# 5. Rule 3. applies until, at some point, intermediate results are not present. From this point onward, all computations will be done. If intermediate results of a later step exist, they will be re-computed.";
	print $fh "# 6. If you are not sure if intermediate results have been computed correctly / completely, please delete the corresponding output file in order to force SANJUAN re-computing these results and all later (down-stream) results.";
	print $fh "";
	print $fh "#######################  General parameters  #######################";
	print $fh "GENOME=hg19      ### ID of genome; run sanjuan -g to see genomes available in this installation";
	print $fh "OUTDIR=          ### Directory for Output, if not given, the current working directory is used as output directory";
	print $fh "COND1=CNTR       ### Label for Condition 1 (eg 'CNT' or 'WT')";
	print $fh "COND2=COND       ### Label for Condition 2 (eg 'KD' or 'OvEx')";
	print $fh "TESTRUN=N        ### values Y, N; if Y, qsub statements are printed but not sent to cluster";
	print $fh "DBLOCATION=      ### the location of the sub-directory db containing predefined exon-exon junctions. Needs to be specified only if db sun-directory is not in main SANJUAN directory.";
	print $fh "NOQSUB=N         ### run SANJUAN without sending jobs to CRG cluster";
	print $fh "NPROCS=1         ### run SANJUAN with this maximal number of parallel processes (only applies to mapping)";	
	print $fh "";
	print $fh "#######################  Data  #######################";
	print $fh "### To specify parameters, un-comment them.";
	print $fh "### Starting with raw FASTQ files: specify the parameter RAWFASTQS to run mapping, and splicing analysis.";
	print $fh "### RAWFASTQS_DIR: directory containing the FASTQ files; These files get assigned to the two groups automatically";
	print $fh "###                and therefore need to fulfill the following naming conventions.";
	print $fh "### Naming convention:";
	print $fh "### 1. GROUP LABELS:     Each file name needs to contain at some point the group label which you defined earlier (COND1 or COND2)";
	print $fh "###                      When choosing lables take care: a label TAI for example is detected (maybe unintentionally) in a file name like 201607_TAIR_R1.fq";
	print $fh "### 2. FILE ENDINGS:     Must be one of these (case-insensitive): fq, fastq, fq.gz, fastq.gz, fq.bz2, fastq.bz2";
	print $fh "### 3. PAIRED-END READS: Must be given in two files, one for reads1 the second for reads2, which comprise a pair of files.";
	print $fh "###                      File names of a pair must be identical except for the last position which must be 1 or 2 indicating reads1, reads2.";
	print $fh "###                      E.g., 201607_exp1_read1.fq and 201607_exp1_read2.fq or SNF1KD_R1.fq.bz2 and SNF1KD_R2.fq.bz2 or control1_1.fq.gz and control1_2.fq.gz";
	print $fh "RAWFASTQS_DIR=          ### folder containing RNAseq; If you start with BAM files, comment this line out with ###.";
	print $fh "";
	print $fh "### If you start with BAM files: specify these parameters BAM1 and BAM2to SKIP mapping, AND ONLY RUN splicing analysis.";
	print $fh "### BAM1 is the merged BAM files for group 1, BAM2 corresponds to group2.";
	print $fh "# BAM1=                  ### Path to condition1 (merged) bam file (contains all mapped reads for this condition)";
	print $fh "# BAM2=                  ### Path to condition2 (merged) bam file (contains all mapped reads for this condition)";
	print $fh "";
	print $fh "#######################  Preprocessing/Mapping parameters  #######################";
	print $fh "ADAPTER=AGATCGGAAGAGC    ### Specify adapter sequence if other than CRG facility default (AGATCGGAAGAGC...).  Leave empty if unknown.";
	print $fh "LIBTYPE=2S               ###library type of RNAseq; values: 1S (Single-end, Stranded), 1U (Single-end, Unstranded), 2U (Paired-end, Unstranded) or 2S (Paired-end, Stranded).Can be omitted if -b B"; #CinS 
	print $fh "PHRED=phred33            ### encoding of base qualities in FASTQ files ASCII+33 -> phred33, ASCII+64 -> phred64 (only used for trimming, see article on the FASTQ file format on Wikipedia)";
	print $fh "TPM=None                 ### None or Basic. Setting to Basic activates the STAR twoPassMode mapping for more sensitive novel junction discovery at the cost of speed.";	#NinS
	print $fh "";
	print $fh "#######################  Splicing Analysis parameters  #######################";
	print $fh "CONF=HC          ### Analysis Stringency Level: 'VHC'-> VeryHighConfidence (DPSI>20%, p-val<0.0001),  'HC'-> HighConfidence (DPSI>15%, p-val<0.001),   'MC'-> MediumConfidence (DPSI>10%, p-val<0.01)"; 
	print $fh "IRM=Y            ### IRM mode: Perform  High Sensitivity Intron Retention Analysis? 'Y'->YES  'N'->NO ";
	print $fh "SUPPJUN=Y        ### Require Supporting Junction Evidence for IntrRet. identification (IRM mode). 'Y'->YES  'N'->NO";
	print $fh "LOWSEQRQMNTS=N   ### Low sequence requirements: set to Y is you are working with RNASeq data not coming from CRG. If set to Y, some stringent tests on RNASeq will be omitted leading to more usable reads.";
	print $fh "RMDUP=N          ### If set to Y: PCR duplicates will be removed\n";
	$\="";
	close($fh);
	exit 0;
}

if(@ARGV==1 && ( $ARGV[0] eq "version" || $ARGV[0] eq "-version" || $ARGV[0] eq "--version") ){
	print "\nSANJUAN $version\n\n";
	exit 0;
}

# check if species is available for splicing analysis
sub is_available_for_splicing{
	my $species=$_[0];
	my $abs_path=$_[1];
	my $check=0;
	# check if we have the genome file
	if(-e "$abs_path/db/genomes/${species}.genome"){$check++;}
	
	# SANJUAN_annotation_files
	if(-e "$abs_path/db/SANJUAN_annotation_files/${species}_Transcript_Junctions.txt" && -e "$abs_path/db/SANJUAN_annotation_files/${species}_Transcripts.bed" && -e "$abs_path/db/SANJUAN_annotation_files/${species}_TxID2Name.txt" ){$check++;}
		
	# gft files are not necessary for splicing analysis
	#if(-e "$abs_path/db/gtfs/${species}.gtf"){$check++;}
		
	# FALSE
	my $ret=0;
	if($check==2){$ret=1;}
return($ret);
}

sub is_available_for_mapping{
	my $species=$_[0];
	my $abs_path=$_[1];
	my $check=0;
	my @txtfiles;
	my @tabfiles;

	# STAR indices
	@txtfiles=<$abs_path/mapping_indexes/$species/*.txt>;
	if(@txtfiles==6){$check++;}
	@tabfiles=<$abs_path/mapping_indexes/$species/*.tab>;
	if(@tabfiles==6){$check++;}
	if(-e "$abs_path/mapping_indexes/$species/Genome"){$check++;}
	if(-e "$abs_path/mapping_indexes/$species/SA"){$check++;}
	if(-e "$abs_path/mapping_indexes/$species/SAindex"){$check++;}	

	# FALSE
	my $ret=0;
	if($check==5){$ret=1;}
return($ret);	
}


if(@ARGV==1 && $ARGV[0] eq "-g"){
	# genome / species shortcuts
	my %all_scs=();
	my %splicing_ok=();
	my $species;
	my %mapping_ok=();

	# available species for splicing analysis 
	foreach my $file (<$abs_path/db/genomes/*.*>){
		if( $file =~ /$abs_path\/db\/genomes\/(.+)\.genome/ ){
			$species=$1;
			if(is_available_for_splicing($species,$abs_path)){$splicing_ok{$species}=1;$all_scs{$species}=1;}
		}
	}
	
	# available species for mapping
	foreach my $dir (<$abs_path/mapping_indexes/*>){
		if( $dir =~ /$abs_path\/mapping_indexes\/(.+)/ ){if(-d "$abs_path/mapping_indexes/$1"){
			$species=$1;
			if(is_available_for_mapping($species,$abs_path)){$mapping_ok{$species}=1;$all_scs{$species}=1;}
		}}
	}

	print "\n   available genomes: ";
	my $str="\n   ID          FOR MAPPING     FOR SPLICING ANALYSIS\n";
	my $str2=$str;
	foreach my $sc (sort {lc $a cmp lc $b} keys %all_scs) {
		$str.= "   $sc";
		for(my $i=3+length($sc);$i<23;$i++){$str.=" ";}
		if(defined($mapping_ok{$sc})){$str.="yes";}else{$str.=" no";}
		for(my $i=0;$i<23;$i++){$str.=" ";}
		if(defined($splicing_ok{$sc})){$str.="yes\n";}else{$str.=" no\n";}
	}

	if($str eq $str2){print "   none\n";}else{print $str;}
	print "\n";

	exit 0;
}


## setting standard values
my $genome="hg19";
my $RNAseq;	#Stranded "S" or unstranded "U" RNAseq experiment. Is set automatically depending on the parameter library_type below. 
my $conf='HC';	#Specify stringency for Differentially Spliced Junctions (VeryHighConfidence -> VHC, HighConfidence -> HC, MediumConfidence -> MC)
my $IRM='N';	#IRM mode: Perform  High Sensitivity Intron Retention Anlaysis. Default 'N' 
my $SuppJun='N';#Require Supporting Junction Evidence for IR identification (IRM mode). Default 'N'
my $output_dir="";
my $phred_code="phred33";#Instructs STAR to interpret fastq quality scores correctly 
my $tpm="None";	#NiS   None or Basic: No two pass mapping by default  
my $library_type="2S";# 1S|1U|2S|2U: Single or Paired end (1|2), Stranded or Unstranded (S|U)
my $adapter="AGATCGGAAGAGC";	#Specify adapter sequence to be trimmed. To find adapter: minion search-adapter -i FASTQFILE.gz (Default adapter sequence for CRG facility, RNAseq )
my $g1_shortname="COND";
my $g2_shortname="CNTR";
# order sample1_read1.fq sample_read2.fq sample2_read1.fq sample2_read2.fq ...
my @g1_files=();# fastq files of group 1
my @g2_files=();# fastq files of group 2
# qsub jobs
my ($bam1,$bam2)=("","","");
my $start_with = "M"; # standard: we start with mapping
my $low_seq_req="N";
my $test_run=0;  # if set to 1, qsub statements will be printed but not sent to cluster
#NRS my $inner_mate_dist=85;
#NRS my $inner_mate_dist_std_dev=25;

# special arguments from parameter file
#  if $map_no_trim=Y -> Go to mapping directly (i.e map using untrimmed fastq files)
#
my $rawinput_dir=".";# input_dir contains all input fastq files 
# this SANJUAN run should not use qsub but run locally
my $run_without_qsub=0;
my $N_processes=1;
my $rmdup=0; #set to 0 to keep PCR duplicates / 1 to discard them

# parameters through arguments
if(@ARGV>1){
	# no standard values for: -f1 -f2 -o 
	for(my $i=0; $i<@ARGV;$i++){
		if($ARGV[$i] eq "-g"){$genome=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-c"){$conf=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-i"){$IRM="Y";}
		if($ARGV[$i] eq "-s"){$SuppJun="Y";}
		if($ARGV[$i] eq "-p"){$phred_code=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-l"){$library_type=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-a"){$adapter=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-g1"){$g1_shortname=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-g2"){$g2_shortname=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-tpm"){$tpm="Basic"}	#NinS
		if($ARGV[$i] eq "-f1"){@g1_files=split(",",$ARGV[($i++)+1]);}
		if($ARGV[$i] eq "-f2"){@g2_files=split(",",$ARGV[($i++)+1]);}
		if($ARGV[$i] eq "-o"){$output_dir=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-b"){$start_with=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-d"){$rawinput_dir=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-lsr"){$low_seq_req="Y";}
		if($ARGV[$i] eq "-t"){$test_run=1;}
		if($ARGV[$i] eq "-db"){$sanjuan_genomic_data_dir=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-noqsub"){$run_without_qsub=1;}
		if($ARGV[$i] eq "-nprocs"){$N_processes=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-rmdup"){$rmdup=1;}
		#NRS if($ARGV[$i] eq "-d"){$inner_mate_dist=$ARGV[($i++)+1];}
		#NRS if($ARGV[$i] eq "-d_dev"){$inner_mate_dist_std_dev=$ARGV[($i++)+1];}
	}
}

# parameters through parameter file
else{	
	# Read and parse parameter file:
	open (my $fh,"<".$ARGV[0]) || die "Cannot open parameter file $ARGV[0] for reading: $!\n";
	while (<$fh>){
		$genome=$1           if $_=~/^\s*GENOME\s*=\s*([\w\/\.\_\-]+)\s*/;	#Species genome: hg-> human, mm-> mouse, dr-> zebrafish
		$rawinput_dir=$1     if $_=~/^\s*RAWFASTQS_DIR\s*=\s*([\w\/\.\_\-\~]+)/;	#Input directory
		#NRS $trimmedinput_dir=$1 if $_=~/^\s*TRIMMEDFASTQS_DIR\s*=([\w\/\.\_\-]+)/;	#Input directory of trimmed fastq files (if given, trimming is skipped)
		$output_dir=$1       if $_=~/^\s*OUTDIR\s*=\s*([\w\/\.\_\-\~]+)/;	#Base directory for Output
		$adapter=$1          if $_=~/^\s*ADAPTER\s*=\s*([ACGTNUacgtnu]+)/;	#Adapter sequence
		$library_type=$1     if $_=~/^\s*LIBTYPE\s*=\s*(1S|1U|2S|2U)/;	#CinS
		$phred_code=$1       if $_=~/^\s*PHRED\s*=\s*(phred33|phred64)/;
		$g1_shortname=$1     if $_=~/^\s*COND1\s*=\s*(\w+)/;
		$g2_shortname=$1     if $_=~/^\s*COND2\s*=\s*(\w+)/;
		$tpm=$1		     if $_=~/^\s*TPM\s*=\s*(Basic|None)/;	#NinS
		$bam1=$1             if $_=~/^\s*BAM1\s*=\s*([\w\/\.\_\-]+)/;	# bam files; if given mapping is skipped
		$bam2=$1             if $_=~/^\s*BAM2\s*=\s*([\w\/\.\_\-]+)/;
		$conf=$1             if $_=~/^\s*CONF\s*=\s*(VHC|HC|MC)/;
		$SuppJun=$1          if $_=~/^\s*SUPPJUN\s*=\s*(Y|N)/;
		$IRM=$1              if $_=~/^\s*IRM\s*=\s*(Y|N)/;
		$low_seq_req=$1      if $_=~/^\s*LOWSEQRQMNTS\s*=\s*(Y|N)/;
		$test_run=1          if $_=~/^\s*TESTRUN\s*=\s*Y/;
		$sanjuan_genomic_data_dir=$1 if $_=~/^\s*DBLOCATION\s*=\s*([\w\/\.\_\-]+)/;
		$run_without_qsub=1  if $_=~/^\s*NOQSUB\s*=\s*Y/;
		$N_processes=$1      if $_=~/^\s*NPROCS\s*=\s*(\d+)/;
		$rmdup=1	     if $_=~/^\s*RMDUP\s*=\s*Y/;
		#NRS $inner_mate_dist=$1  if $_=~/^\s*INNER_MATE_DIST\s*=(\d+)/;
		#NRS $inner_mate_dist_std_dev=$1  if $_=~/^\s*DIST_STD_DEV\s*=(\d+)/;
	}
	close($fh);
}


## parameter checks
my ($OK_params_preprocess,$OK_params_main)=(1,1);
my ($warnings_preprocess,$warnings_main,$tmp_str)=("","",);
unless($genome){$tmp_str="Parameter GENOME/-g not defined.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess.=$tmp_str;$warnings_main.=$tmp_str;}
unless($conf =~ /VHC|HC|MC/){$tmp_str="Parameter CONF/-c not or wrongly defined. Should take values VHC, HC, MC.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
unless($IRM =~ /Y|N/){$tmp_str="Parameter IRM/-i not or wrongly defined. Should take values Y or N.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
unless($SuppJun =~ /Y|N/){$tmp_str="Parameter SUPPJUN/-s not or wrongly defined. Should take values Y or N.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
unless($phred_code =~ /phred33|phred64/){$tmp_str="Parameter PHRED/-p not or wrongly defined. Should take values phred33 or phred64.\n";$OK_params_preprocess=0;$warnings_preprocess.=$tmp_str;}
unless($library_type =~ /1S|1U|2S|2U/)  {$tmp_str="Parameter LIBTYPE/-l not defined. Should take values 1S,1U, 2S or 2U\n"; $OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess.=$tmp_str;$warnings_main.=$tmp_str;}
unless($tpm =~ /Basic|None/)  {$tmp_str="Parameter TPM/-tpm not defined. Should take values Basic or None\n"; $OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess.=$tmp_str;$warnings_main.=$tmp_str;}
unless($adapter =~ /[ACGTNUacgtnu]+/){$tmp_str="Parameter ADAPTER/-a not defined. Should be a sequence composed of any letter of ACGTNUacgtnu.\n";$OK_params_preprocess=0;$warnings_preprocess.=$tmp_str;}
unless($g1_shortname =~ /\w+/){$tmp_str="Parameter COND1/-g1 not defined. Should be a short word composed of a-z, A-Z, 0-9 and \_.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess.=$tmp_str;$warnings_main.=$tmp_str;}
unless($g2_shortname =~ /\w+/){$tmp_str="Parameter COND2/-g2 not defined. Should be a short word composed of a-z, A-Z, 0-9 and \_.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess.=$tmp_str;$warnings_main.=$tmp_str;}
if($g1_shortname eq $g2_shortname){$tmp_str="Parameter COND1/-g1 and COND2/-g2 are identical but must be different.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess.=$tmp_str;$warnings_main.=$tmp_str;}
unless($low_seq_req =~ /Y|N/){$tmp_str="Parameter LOWSEQRQMNTS/-r not or wrongly defined. Should take values Y or N.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
#NRS unless($inner_mate_dist =~ /(\d+)/){$tmp_str="Parameter INNER_MATE_DIST/-d not or wrongly defined. Should take one integer value.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
#NRS unless($inner_mate_dist_std_dev =~ /(\d+)/){$tmp_str="Parameter DIST_STD_DEV/-d_dev not or wrongly defined. Should take one integer value.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}


if($output_dir eq "" || $output_dir eq "."){$output_dir=cwd();}
# remove trailing / if exists
$rawinput_dir=~ s/\/$//;
$output_dir=~ s/\/$//;


if(@g1_files>0){
	unless($start_with eq "M" || $start_with eq "S" ){$tmp_str="Parameter -b not defined. Should be set to M or S to start with mapping, or directly with splicing analysis\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess.=$tmp_str;$warnings_main.=$tmp_str;} #CinS
	# add data input dir if necessary
	unless($rawinput_dir eq "."){
		for(my $i=0;$i<@g1_files;$i++){$g1_files[$i]=$rawinput_dir."/".$g1_files[$i];}
		for(my $i=0;$i<@g2_files;$i++){$g2_files[$i]=$rawinput_dir."/".$g2_files[$i];}
	}
	
	if($start_with eq "S"){
		my $check=0;
		my $tmp_msg="";
		if(@g1_files>1){$tmp_msg.="More than one input file for group 1 given. If parameter -b is set to S you want to start directly with the splicing analysis and you should define only one BAM input file for group 1.\n";$check=1;}
		if(@g2_files>1){$tmp_msg.="More than one input file for group 2 given. If parameter -b is set to S you want to start directly with the splicing analysis and you should define only one BAM input file for group 2.\n";$check=1;}
		if($check){die $tmp_msg;}
		
		$bam1=$g1_files[0];
		$bam2=$g2_files[0];
	}
}else{
	unless(($bam1 && $bam2) || $rawinput_dir){
		die "No input is specified. Specify either BAM1/BAM2 or RAWFASTQS_DIR\n";
	}

	if($bam1 && $bam2){
		$start_with="S"; # directly go to splicing analysis
	}else{# start with mapping; standard
		$start_with="M";
		
		# generate arrays with fastq files for group1 and group2 if given through parameter file
		if($start_with eq "M" && @g1_files==0 && @g2_files==0){	#CinS
			my $input_dir = $rawinput_dir; #CinS
			# file names should be like *_[g1shortname|g2shortname]_*_[r|read|R][1|2].[fq|fastq].[gz]
			# .gz is optional
			# files get sorted into the two groups according to [g1shortname|g2shortname] and according to read1/2 according to [r|read|R][1|2]  
			foreach (<$input_dir/*>){
				# file name with full path
				my $file=$_;
				next unless (-f $file);

				if($file =~ /$g1_shortname.*\.(fq$|fastq$|fq\.gz$|fastq\.gz$|fq\.bz2$|fastq\.bz2$)/i){
					push(@g1_files,$file);
					next;
				}
				if($file =~ /$g2_shortname.*\.(fq$|fastq$|fq\.gz$|fastq\.gz$|fq\.bz2$|fastq\.bz2$)/i){
					push(@g2_files,$file);
					next;
				}
			}
			# files with path
			@g1_files=sort(@g1_files);
			@g2_files=sort(@g2_files);
		}		
	}
}

# stop if something is wrong with parameters
if($start_with ne "S" && !$OK_params_preprocess){die $warnings_preprocess;}
if(!$OK_params_main){die $warnings_main;}

# is set automatically depending on library type
$RNAseq=($library_type=~/S/i)? "S":"U";	#CinS



# @g1_files and @g2_files might contain two different kinds of files depending on start_with 
# start_with=T(rimming) -> they contain raw fastq files
# start_with=M(apping) -> they contain trimmed fastq files
# start_with=S (splicing analysis) -> they contain for each group one bam file
#
# file checks: do files exist and can they be opened


if($start_with eq "S"){
	open(my $fh,"<".$bam1) or die "Cannot open BAM file $bam1: $!\n";close($fh);
	open($fh,"<".$bam2) or die "Cannot open BAM file $bam2: $!\n";close($fh);
}else{
	if(@g1_files <1 || ($library_type=~/2/ && scalar(@g1_files) % 2 != 0)){die "Wrong number of FASTQ input files for group 1. Number must be non-zero and even for paired-end read data.\n";}

	if(@g2_files <1 || ($library_type=~/2/ && scalar(@g2_files) % 2 != 0)){die "Wrong number of FASTQ input files for group 2. Number must be non-zero and even for paired-end read data.\n";}
	foreach my $fname (@g1_files,@g2_files){
		# do files exsist and can be opened?
		open(my $fh,"<".$fname) or die "Cannot open FASTQ file $fname: $!\n";close($fh);
		# check file endings 
		unless( $fname =~ /(\.fq\.bz2$|fastq\.bz2$|\.fastq\.gz$|\.fq\.gz$|\.fq$|\.fastq$)/i ){die "FASTQ input files should have endings fastq, fq, fastq.gz, fq.gz, fastq.bz2 or fq.bz2 but file $fname doesn't have\n";}
	}
}



# tophat files:
# transcriptome index, gene annotation, bowtie index
#NRS my ($tophat_tr_index,$tophat_bowtie_index)=($abs_path."/mapping_indexes/$genome/${genome}_transcriptome",$abs_path."/mapping_indexes/$genome/$genome");
my $STAR_index=$abs_path."/mapping_indexes/$genome/"; #NiS

if($start_with ne "S"){
	# user wants to do mapping and splicing analysis
	if(!is_available_for_mapping($genome,$abs_path)){die "Mapping for species / genome $genome is not possible\nas this species / genome is not installed for mapping.\n Run sanjuan -g to see installed species / genomes, \nsee README on GitHub.com webpage of SANJUAN on how to install further species / genomes\n";}
}else{
	# user wants to do splicing analysis only
	if(!is_available_for_splicing($genome,$abs_path)){die "Splicing analysis for species / genome $genome is not possible\nas this species / genome is not installed for splicing analysis.\n Run sanjuan -g to see installed species / genomes, \nsee README on GitHub.com webpage of SANJUAN on how to install further species / genomes\n";}
}


#########
#########
my $call;
print "*************\nSANJUAN\n*************\n\n";

my $suffix="";
if($test_run){$suffix.="-t ";}
if($run_without_qsub){$suffix.="-nqsub ";}

print "Call:\nsanjuan -g $genome -c $conf -nproc $N_processes -i $IRM -s $SuppJun -l $library_type -a $adapter -tpm $tpm -g1 $g1_shortname -f1 ".join(",",@g1_files)." -g2 $g2_shortname -f2 ".join(",",@g2_files)." -o $output_dir -b $start_with -r $low_seq_req $suffix\n\n\n";

unless (-d $output_dir){print `mkdir -p $output_dir`;}
# here go all output and error messages
unless (-d "$output_dir/log_files"){print `mkdir -p $output_dir/log_files`;}
my $ret="Job ids:";
if($start_with eq "M"){
	# mapping
	print "\n\n*************\nMapping\n*************\n\n";
	$call="perl $sanjuan_dir/preProcess_and_Map.pl $output_dir $genome $adapter $library_type $tpm $test_run $run_without_qsub $sanjuan_dir $STAR_index $N_processes $phred_code -g1 $g1_shortname @g1_files -g2 $g2_shortname @g2_files"; #CinS
	$ret=`$call`;
	print $ret."\n";
}else{
	print "\n\n\n\nMapping\n#####################\n\nSkipped\n\n\n\n\nMerging BAM files\n#####################\n\nSkipped\n\n";
}

my $job_ids=-1;
my @fs=split("\n",$ret);
if(@fs>0){
	chomp($fs[@fs-1]);
	$job_ids=$1 if $fs[@fs-1]=~/Job ids:(.+)/;
}

print "job ids from preProcess_and_Map: $job_ids\n\n\n";


my $merged_bam_file_1 = ($start_with eq "S")? $bam1 : $output_dir . '/MAPPING/' . $g1_shortname ."_Aligned.sortedByCoord.out.merged.bam";
my $merged_bam_file_2 = ($start_with eq "S")? $bam2 : $output_dir . '/MAPPING/' . $g2_shortname ."_Aligned.sortedByCoord.out.merged.bam";

print "merged_bam_file_1=$merged_bam_file_1\nmerged_bam_file_2=$merged_bam_file_2\n\n";



print "\n\n*************\nSplicing Analysis\n*************\n\n";
$call="perl $sanjuan_dir/SANJUAN_wrapper.pl $genome $library_type $conf $IRM $SuppJun $g1_shortname $g2_shortname $merged_bam_file_1 $merged_bam_file_2 $output_dir $job_ids $low_seq_req $test_run $run_without_qsub $N_processes $rmdup $sanjuan_dir $sanjuan_perllib $sanjuan_genomic_data_dir";
system($call);

exit(0);
