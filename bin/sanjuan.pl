#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

### Verion
my $version="1.0";

#### main paths parameters
# is the full path to the script sanjuan.pl including sanjuan.pl
# even, if sanjuan was called through a link 
my $abs_path=abs_path($0);
# remove /bin/sanjuan.pl to obtain root directory of SANJUAN installation
$abs_path =~ s/\/bin\/sanjuan\.pl//;
my $sanjuan_dir=$abs_path."/lib";   			# should contain all scripts used by SANJUAN
my $sanjuan_perllib=$abs_path."/perllib";		# should contain Perl modules Text and Statistics XXX check if we need both!
my $sanjuan_genomic_data_dir=$abs_path."/db";	# should contain folders genomes and annotation_files

my $user="unpecified"; # switch between Andre and Pan to overwrite path variables
if($user eq "Andre"){	
	$sanjuan_dir="/users/mirimia/agohr/crg/projects/2015_sanjuan/git/SANJUAN/lib";
	$sanjuan_perllib="/users/mirimia/agohr/crg/projects/2015_sanjuan/git/SANJUAN/perllib";
	$sanjuan_genomic_data_dir="/users/mirimia/agohr/crg/projects/2015_sanjuan/git/SANJUAN/db";
}
if($user eq "Pan"){
	$sanjuan_dir="";
	$sanjuan_perllib=""; 
	$sanjuan_genomic_data_dir="";  
}

# Are all important parts of SANJUAN in place?
# check if all sanjuan files can be accessed
sub check_file_access{open(my $fh,"<".$_[0]) or die "Cannot open file $_[0] which is essential for SANJUAN. You might go through the installation process again to solve this problem.\n";close($fh);}
my @SANJUAN_files=("annotate_Diff_Used_Introns.pl","annotate_Diff_Used_Junctions.pl","calc_INTRON_retention.pl","calc_JUNCT_efficiency.pl","get_juncts.pl","job1.pl","job2.pl","ln_s_wrapper.sh","merge_junctions.pl","preProcess_and_Map.pl","SANJUAN_wrapper.pl","sort_wrapper.sh","tophat_junctions2IntronSegments.pl");
foreach (@SANJUAN_files){check_file_access($sanjuan_dir."/".$_);}

# check if sub-directories genomes and annotation_files exist in db sub-directory
if(!-d $sanjuan_genomic_data_dir."/genomes" || !-d $sanjuan_genomic_data_dir."/annotation_files"){die "Sud-directories genomes and/or annotation_files cannot be found in $sanjuan_genomic_data_dir. Try to specify parameter DBLOCATION / -db to specify their location\n";}

# check if sub-directories perllib is in place
if(!-d $sanjuan_perllib."/Text" || !-d $sanjuan_perllib."/Statistics"){die "Sud-directories / perl modules Text and/or Statistics cannot be found in $sanjuan_perllib. You might go through the installation process again to solve this problem.\n";}
#############################################################
#############################################################


sub print_help{
	print "\nSANJUAN $version -- *S*plicing *AN*alysis & *JU*nction *AN*notation\n===========\n\n";
	print "Parameter definition by a dedicated parameter file or by arguments of the program call.\n\n";
	print "Call with parameter file: sanjuan <parameter-file>\n";
	print "To create a standard parameter file SANJUAN_parameters.txt with explanations run: sanjuan -exampleF\n";
	print "\n";
	print "Call with arguments (<standard value already set>):\n";
	print "sanjuan -g <hg> -g1 <grp1> -f1 -g2 <grp2> -f2 -o <.> -p <phred33> -l <fr-firststrand> -a <AGATCGGAAGAGC> -b <T> -c <HC> -i -s -lsr -t -db <db sub-directory of SANJUAN directory>\n\n";
	print "\t-g:    genome / species; values: hg -> human, mm-> mouse, dr-> zebrafish\n";
	print "\t-g1:   short name for group 1\n";
	print "\t-f1:   input files for group 1; Depending on value of argument -b, -f1 defines different input files.\n";
	print "\t\t if -b is T or M: pairs of FASTQ/FASTQ.GZ files describing paired-end RNAseq data,\n";
	print "\t\t\t comma separated list without white spaces following this order (file names don't matter): run1_read1.fastq,run1_read2.fastq,run2_read1.fastq,run2_read2.fastq,..\n";
	print "\t\t if -b set to B: exactly one BAM file containing all mapped reads for group 1\n";
	print "\t-g2:   like argument -g1 but for group 2\n";
	print "\t-f2:   like argument -f1 but for group 2\n";
	print "\t-o:    output directory; If omitted current working directory is taken.\n";
	print "\t-p:    encoding of fastq qualities; values phred33 -> ASCII+33, phred64 -> ASCII+64; Can be omitted if -b B.\n";
	print "\t\t For details have a look at section Encoding of Wikipedia article on FASTQ https://en.wikipedia.org/wiki/Fastq\n";	
	print "\t-l:    library type of RNAseq samples; values fr-unstranded, fr-firststrand, fr-secondstrand (fr-firststrand is standard for CRG samples).  Can be omitted if -b B\n";
	print "\t\t Set to fr-unstranded if you are in doubt about the library type; though this parameter is critical and it is recommended\n";
	print "\t\t setting it correctly according to the RNAseq data you are working with.\n";
	print "\t\t For details have a look at section on library type of TopHat online manual https://ccb.jhu.edu/software/tophat/manual.shtml.\n";
	print "\t-a:    adapter used in RNAseq measurement; CRG standard is AGATCGGAAGAGC. Can be omitted if -b B.\n"; 
	print "\t\t To identify adapter you could try minion search-adapter -i FASTQFILE.gz)\n";
	print "\t-b:    Starting point; values T -> start with trimming, M -> start with mapping SKIPPING trimming, B- > start with splicing analysis SKIPPING trimming and mapping\n";	
	print "\t-c:    threshold on reported differentially spliced junctions; values VHC -> very high confidence (DPSI>20%, p-val<0.0001),\n";
	print "\t\t HC -> high confidence (DPSI>15%, p-val<0.001), MC -> medium confidence (DPSI>10%, p-val<0.01)\n";
	print "\t-i:    If -i is given, high sensitivity intron retention anlaysis (IRM mode) will be done.\n";
	print "\t-s:    If -s is given, supporting junction evidence for IR identification (for IRM mode -i) will be required; not required if -s omitted\n";
	print "\t-lsr:  If -lsr (low sequence requirements) is given, reads will be filter-out less strictly.\n";
	print "\t-t:    If -t (test run) is given, qsub statements will be printed but not sent to cluster.\n";
	print "\t-db: Full path to the directory db where SANJUAN will find pre-defined exon-exon junctions, genomes, and annotations.\n";
	print "\t\t Has to be set only if this directory is not under the SANJUAN installation directory.\n";
	print "\nFor printing a full example SANJUAN call: sanjuan -exampleC\n\n";
	print "Example call for human RNAseq data from CRG:\n";
	print "\tsanjuan -g1 ko -g2 cntr -f1 run1_1.fastq,run1_2.fastq -f2 run2_1.fastq,run2_2.fastq -c HC -i -s\n\n";
	print "Contact: Panagiotis Papasaikas\n\n";
}

if(@ARGV==0 || $ARGV[0] eq "--help" || $ARGV[0] eq "-help" || $ARGV[0] eq "help" || $ARGV[0] eq "?"){
	print_help;
	exit(0);
}

if(@ARGV==1 && $ARGV[0] eq "-exampleC"){
	print "\nsanjuan -g hg -g1 grp1 -f1 <file1,file2,..> -g2 grp2 -f2 <file1,file2,..> -p phred33 -l fr-firststrand -a AGATCGGAAGAGC -b T -c HC -i -s -lsr -t -db\n\n";
	exit 0;
}

if(@ARGV==1 && $ARGV[0] eq "-exampleF"){
	if(-f "SANJUAN_parameters.txt"){print "File SANJUAN_parameters.txt already exists. Will not write example parameter file.\n";
		exit 0;
	}
	open(my $fh,">"."SANJUAN_parameters.txt") or die "Cannot open file SANJUAN_parameters.txt for writing";

	$\="\n";	
	print $fh "#######################  General notes on behavior of pipeline  #######################";
	print $fh "# 	1. SANJUAN has three starting points from which you can start the analysis depending on the data you are working with.";
	print $fh "# 1.a raw FASTQ files: specify the parameter RAWFASTQS in section Data to run trimming, mapping, and splicing analysis.";
	print $fh "# 1.b trimmed FASTQ files: specify the parameter TRIMMEDFASTQS in section Data to SKIP trimming, AND ONLY RUN mapping and splicing analysis.";
	print $fh "# 1.c BAM files: specify the parameters BAM1 and BAM2 in section Data to SKIP trimming and mapping, AND ONLY RUN splicing analysis.";
	print $fh "# 2. SANJUAN will always start at the lates possible starting point. E.g., if you specify BAM1 and BAM2, SANJUAN will start directly with the splicing analysis regardless of whether you have specified TRIMMEDFASTQS or RAWFASTQS or none of both.";
	print $fh "#";
	print $fh "# 3. To avoid un-necessary computations, SANJUAN applies the following rules to decide which computations are skipped further.";
	print $fh "# 4. If intermediate output files / results of preprocessing or splicing analysis are present, they will be used and not re-computed.";
	print $fh "# 5. Rule 3. applies until, at some point, intermediate results are not present. From this point onward, all computations will be done. If intermediate results of a later step exist, they will be re-computed.";
	print $fh "# 6. If you are not sure if intermediate results have been computed correctly / completely, please delete the corresponding output file in order to make SANJUAN re-computing these results and all later (down-stream) results.";
	print $fh "";
	print $fh "";
	print $fh "#######################  General parameters  #######################";
	print $fh "GENOME=hg		### Specify Organism: hg-> human, mm-> mouse, dr-> zebrafish";
	print $fh "OUTDIR=/users/jvalcarcel/ppapasaikas/SOPHIE/SF3B1_Ast_Data/		### Directory for Output, if not given, the current working directory is used as output directory";
	print $fh "COND1=CNT		### Label for Condition 1 (eg 'CNT' or 'WT')";
	print $fh "COND2=KD			### Label for Condition 2 (eg 'KD' or 'OvEx')";
	print $fh "TESTRUN=N		### values Y, N; if Y, qsub statements are printed but not sent to cluster";
	print $fh "DBLOCATION=		### the location of the sub-directory db containing predefined exon-exon junctions. Needs to be specified only if db sun-directory is not in main SANJUAN directory.\n";
	print $fh "";
	print $fh "";
	print $fh "#######################  Data  #######################";
	print $fh "### To specify parameters, uncomment them.";
	print $fh "### You start with raw FASTQ files: specify this parameter RAWFASTQS to run trimming, mapping, and splicing analysis.";
	print $fh "### Directory containing the FASTQ files";
	print $fh "### Naming convention:\n";
	print $fh "### 1. file names have to contain at some point the label COND1 or COND2 according to which they get assigned to conditions\n";
	print $fh "### 2. endings should be (r|read|R)(1|2).(fq|fastq|fq.gz|fastq.gz)\n";
	print $fh "### 3. the only difference between the read-pair files should be read1 vs read2 or r1 vs r2 or R1 vs R2\n";
	print $fh "### examples for COND1=cntr: test_cntr_somethingmore_r1.fq or siRNA_cntr_something_R2.fq.gz or  dec_cntr_read2.fastq.gz\n";
	print $fh "RAWFASTQS_DIR=/users/jvalcarcel/ppapasaikas/SOPHIE/SF3B1_Ast_Data/";
	print $fh "";
	print $fh "### You start with trimmed FASTQ files: specify this parameter TRIMMEDFASTQS to SKIP trimming, AND ONLY RUN mapping and splicing analysis.";
	print $fh "### TRIMMEDFASTQS_DIR=/users/jvalcarcel/ppapasaikas/SOPHIE/SF3B1_Ast_Data/";
	print $fh "";
	print $fh "### You start with BAM files: specify these parameters BAM1 and BAM2to SKIP trimming and mapping, AND ONLY RUN splicing analysis.";
	print $fh "### BAM1 is the joined BAM files for group 1, BAM2 corresponds to group2.";
	print $fh "# BAM1=/users/jvalcarcel/ppapasaikas/SOPHIE/SF3B1_Ast_Data/TOPHAT_NTsi/accepted_hits.bam		### Path to condition1 (merged) bam file";
	print $fh "# BAM2=/users/jvalcarcel/ppapasaikas/SOPHIE/SF3B1_Ast_Data/TOPHAT_SF3B1/accepted_hits.bam		### Path to condition2 (merged) bam file";
	print $fh "";
	print $fh "";
	print $fh "#######################  Preprocessing and Mapping parameters  #######################";
	print $fh "ADAPTER=AGATCGGAAGAGC	### Specify adapter sequence if other than CRG facility default (AGATCGGAAGAGC...).  Leave empty if unknown.";
	print $fh "LIBTYPE=fr-firststrand	### Specify sequencing library type (fr-firststrand, fr-secondstrand, fr-unstranded). See tophat online manual for details. fr-firstrsand default for CRG facility. Leave empty if unknown.";
	print $fh "PHRED=phred33		### encoding of fastq qualities ASCII+33 -> phred33, ASCII+64 -> phred64 (only used for trimming, see details on wikipedia for fastq file format)";
	print $fh "";
	print $fh "";
	print $fh "#######################  Splicing Analysis parameters  #######################";
	print $fh "CONF=HC			### Analysis Stingency Level: 'VHC'-> VeryHighConfidence (DPSI>20%, p-val<0.0001),  'HC'-> HighConfidence (DPSI>15%, p-val<0.001),   'MC'-> MediumConfidence (DPSI>10%, p-val<0.01)"; 
	print $fh "IRM=Y			### IRM mode: Perform  High Sensitivity Intron Retention Analysis? 'Y'->YES  'N'->NO ";
	print $fh "SUPPJUN=Y		### Require Supporting Junction Evidence for IntrRet. identification (IRM mode). 'Y'->YES  'N'->NO";
	print $fh "LOWSEQRQMNTS=N		### Low sequence requirements: set to Y is you are working with RNASeq data not coming from CRG. If set to Y, some stringent tests on RNASeq will be omitted leading to more usable reads.";
	$\="";
	
	close($fh);
	exit 0;
}


## setting standard values
my $genome="hg";# hg -> human, mm-> mouse, dr-> zebrafish
my $RNAseq="S";	#Stranded "S" or unstranded "U" RNAseq experiment. Is set automatically depending on the parameter library_type (fr-first/secondstrand -> "S", fr-unstranded -> "U"). 
my $conf='HC';	#Specify stringency for Differentially Spliced Junctions (VeryHighConfidence -> VHC, HighConfidence -> HC, MediumConfidence -> MC)
my $IRM='N';	#IRM mode: Perform  High Sensitivity Intron Retention Anlaysis. Default 'N' 
my $SuppJun='N';#Require Supporting Junction Evidence for IR identification (IRM mode). Default 'N'
my $output_dir="";
my $phred_code="phred33";#Instructs Cutadapt to use ASCII+33 / ASCII+64 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding) / (Illumina 1.5 encoding) for quality trimming.
my $library_type="fr-firststrand";# fr-unstranded, fr-firststrand, fr-secondstrand, http://onetipperday.blogspot.com.es/2012/07/how-to-tell-which-library-type-to-use.html
my $adapter="AGATCGGAAGAGC";	#Specify adapter sequence to be trimmed. To find adapter: minion search-adapter -i FASTQFILE.gz (Default adapter sequence for CRG facility, RNAseq )
my $g1_shortname="grp1";
my $g2_shortname="grp2";
# order sample1_read1.fq sample_read2.fq sample2_read1.fq sample2_read2.fq ...
my @g1_files=();# fastq files of group 1
my @g2_files=();# fastq files of group 2
# qsub jobs
my ($bam1,$bam2)=("","","");
my $start_with = "T"; # standard: we start with trimming
my $low_seq_req="N";
my $test_run=0;  # if set to 1, qsub statements will be printed but not sent to cluster

# special arguments from parameter file
#  if $map_no_trim=Y -> Go to mapping directly (i.e map using untrimmed fastq files)
#
my ($rawinput_dir,$trimmedinput_dir)=("","");# input_dir contains all input fastq files 

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
		if($ARGV[$i] eq "-f1"){@g1_files=split(",",$ARGV[($i++)+1]);}
		if($ARGV[$i] eq "-f2"){@g2_files=split(",",$ARGV[($i++)+1]);}
		if($ARGV[$i] eq "-o"){$output_dir=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-b"){$start_with=$ARGV[($i++)+1];}
		if($ARGV[$i] eq "-lsr"){$low_seq_req="Y";}
		if($ARGV[$i] eq "-t"){$test_run=1;}
		if($ARGV[$i] eq "-db"){$sanjuan_genomic_data_dir=$ARGV[($i++)+1];}
	}
}

# parameters through parameter file
else{	
	# Read and parse parameter file:
	open (my $fh,"<".$ARGV[0]) || die "Cannot open parameter file $ARGV[0] for reading: $!\n";
	while (<$fh>){
		$genome=$1 if $_=~/^\s*GENOME\s*=(hg|mm|dr)/;	#Species genome: hg-> human, mm-> mouse, dr-> zebrafish
		$rawinput_dir=$1 if $_=~/^\s*RAWFASTQS_DIR\s*=([\w\/\.\_\-]+)/;	#Input directory
		$trimmedinput_dir=$1 if $_=~/^\s*TRIMMEDFASTQS_DIR\s*=([\w\/\.\_\-]+)/;	#Input directory of trimmed fastq files (if given, trimming is skipped)
		$output_dir=$1 if $_=~/^\s*OUTDIR\s*=([\w\/\.\_\-]+)/;	#Base directory for Output
		$adapter=$1 if $_=~/^\s*ADAPTER\s*=([ACGTNUacgtnu]+)/;	#Adapter sequence
		$library_type=$1 if $_=~/^\s*LIBTYPE\s*=(fr\-firststrand|fr\-secondstrand|fr\-unstranded)/;
		$phred_code=$1 if  $_=~/^\s*PHRED\s*=(phred33|phred64)/;
		$g1_shortname=$1 if  $_=~/^\s*COND1\s*=(\w+)/;
		$g2_shortname=$1 if  $_=~/^\s*COND2\s*=(\w+)/;
		$bam1=$1 if $_=~/^\s*INFDIR\s*=([\w\/\.\_\-]+)/;	# bam files; if given trimming and mapping are skipped
		$bam2=$1 if $_=~/^\s*INFDIR\s*=([\w\/\.\_\-]+)/;
		$conf=$1 if $_=~/^\s*CONF\s*=(VHC|HC|MC)/;
		$SuppJun=$1 if $_=~/^\s*SUPPJUN\s*=(Y|N)/;
		$IRM=$1 if $_=~/^\s*IRM\s*=(Y|N)/;
		$low_seq_req=$1 if $_=~/^\s*LOWSEQRQMNTS\s*=(Y|N)/;
		$test_run=1 if $_=~/^\s*TESTRUN\s*=Y/;
		$sanjuan_genomic_data_dir=$1 if $_=~/^\s*DBLOCATION\s*=([\w\/\.\_\-]+)/;
	}
	close($fh);
}

## parameter checks
my ($OK_params_preprocess,$OK_params_main)=(1,1);
my ($warnings_preprocess,$warnings_main,$tmp_str)=("","",);
unless($genome =~ /hg|mm|dr/){$tmp_str="Parameter GENOME/-g not or wrongly defined. Should take values hg, mm, dr.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess-=$tmp_str;$warnings_main.=$tmp_str;}
unless($conf =~ /VHC|HC|MC/){$tmp_str="Parameter CONF/-c not or wrongly defined. Should take values VHC, HC, MC.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
unless($IRM =~ /Y|N/){$tmp_str="Parameter IRM/-i not or wrongly defined. Should take values Y or N.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
unless($SuppJun =~ /Y|N/){$tmp_str="Parameter SUPPJUN/-s not or wrongly defined. Should take values Y or N.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}
unless($phred_code =~ /phred33|phred64/){$tmp_str="Parameter PHRED/-p not or wrongly defined. Should take values phred33 or phred64.\n";$OK_params_preprocess=0;$warnings_preprocess-=$tmp_str;}
unless($library_type =~ /fr\-firststrand|fr\-secondstrand|fr\-unstranded/){$tmp_str="Parameter LIBTYPE/-l not defined. Should take values fr-firststrand, fr-secondstrand, fr-unstranded.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess-=$tmp_str;$warnings_main.=$tmp_str;}
unless($adapter =~ /[ACGTNUacgtnu]+/){$tmp_str="Parameter ADAPTER/-a not defined. Should be a sequence composed of any letter of ACGTNUacgtnu.\n";$OK_params_preprocess=0;$warnings_preprocess-=$tmp_str;}
unless($g1_shortname =~ /\w+/){$tmp_str="Parameter COND1/-g1 not defined. Should be a short word composed of a-z, A-Z, 0-9 and \_.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess-=$tmp_str;$warnings_main.=$tmp_str;}
unless($g2_shortname =~ /\w+/){$tmp_str="Parameter COND2/-g2 not defined. Should be a short word composed of a-z, A-Z, 0-9 and \_.\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess-=$tmp_str;$warnings_main.=$tmp_str;}
unless($low_seq_req =~ /Y|N/){$tmp_str="Parameter LOWSEQRQMNTS/-r not or wrongly defined. Should take values Y or N.\n";$OK_params_main=0;$warnings_main.=$tmp_str;}

if(@g1_files>0){
	unless($start_with eq "M" || $start_with eq "T" || $start_with eq "B" ){$tmp_str="Parameter -b not defined. Should be set to T or M or B to start with trimming, mapping, or directly with splicing analysis\n";$OK_params_preprocess=0;$OK_params_main=0;$warnings_preprocess-=$tmp_str;$warnings_main.=$tmp_str;}
	if($start_with eq "B"){
		my $check=0;
		my $tmp_msg="";
		if(@g1_files>1){$tmp_msg.="More than one input file for group 1 given. If parameter -b is set to B you want to start directly with the splicing analysis and you should define only one BAM input file for group 1.\n";$check=1;}
		if(@g2_files>1){$tmp_msg.="More than one input file for group 2 given. If parameter -b is set to B you want to start directly with the splicing analysis and you should define only one BAM input file for group 2.\n";$check=1;}
		if($check){die $tmp_msg;}
		
		$bam1=$g1_files[0];
		$bam2=$g2_files[0];
	}
}else{
	unless(($bam1 && $bam2) || $rawinput_dir || $trimmedinput_dir){
		die "No input is specified. Specify either BAM1/BAM2 or RAWFASTQS_DIR or TRIMMEDFASTQS_DIR\n";
	}

	if($bam1 && $bam2){$start_with="B"; # directly to go splicing analysis
	}elsif($trimmedinput_dir){$start_with="M" # start with mapping, skip trimming
	}else{$start_with="T"}  # start with trimming; standard
}

# stop if something is wrong with parameters
if($start_with ne "B" && !$OK_params_preprocess){die $warnings_preprocess;}
if(!$OK_params_main){die $warnings_main;}

# is set automatically depending on library type
if($library_type eq "fr-unstranded"){$RNAseq="U";}else{$RNAseq="S";}

if($output_dir eq "" || $output_dir eq "."){$output_dir=cwd();}
# remove trailing / if exists
$rawinput_dir=~ s/\/$//;
$trimmedinput_dir=~ s/\/$//;
$output_dir=~ s/\/$//;

# generate arrays with fastq files for group1 and group2 if given through parameter file
if($start_with ne "B" && @g1_files==0 && @g2_files==0){
	my $input_dir = ($start_with eq "M")? $trimmedinput_dir : $rawinput_dir;
	# file names should be like *_[g1shortname|g2shortname]_*_[r|read|R][1|2].[fq|fastq].[gz]
	# .gz is optional
	# files get sorted into the two groups according to [g1shortname|g2shortname] and according to read1/2 according to [r|read|R][1|2]  
	foreach (<$input_dir/*>){
		# file name with full path
		my $file=$_;
		next unless (-f $file);
		
		if($file =~ /$g1_shortname.*(r|read|R)(1|2)\.(fq$|fastq$|fq\.gz$|fastq\.gz$)/){
			push(@g1_files,$file);
			next;
		}
		if($file =~ /$g2_shortname.*(r|read|R)(1|2)\.(fq$|fastq$|fq\.gz$|fastq\.gz$)/){
			push(@g2_files,$file);
			next;
		}
	}
	# files with path
	@g1_files=sort(@g1_files);
	@g2_files=sort(@g2_files);
}

# @g1_files and @g2_files might contain two different kinds of files depending on start_with 
# start_with=T(rimming) -> they contain raw fastq files
# start_with=M(apping) -> they contain trimmed fastq files
# start_with=B (splicing analysis) -> they contain for each group one bam file
#
# file checks: do files exist and can they be opened
if($start_with eq "B"){
	open(my $fh,"<".$bam1) or die "Cannot open BAM file $bam1: $!\n";close($fh);
	open($fh,"<".$bam2) or die "Cannot open BAM file $bam2: $!\n";close($fh);
}else{
	if(@g1_files <2 || scalar(@g1_files) % 2 != 0){die "Wrong number of FASTQ input files for group 1 given. Number must be even and at least 2.\n";}
	if(@g2_files <2 || scalar(@g2_files) % 2 != 0){die "Wrong number of FASTQ input files for group 2 given. Number must be even and at least 2.\n";}	
	foreach my $fname (@g1_files,@g2_files){
		# do files exsist and can be opened?
		open(my $fh,"<".$fname) or die "Cannot open FASTQ file $fname: $!\n";close($fh);
		# check file endings 
		unless( $fname =~ /(\.fastq\.gz$|\.fq\.gz$|\.fq$|\.fastq$)/ ){die "FASTQ input files should have endings fastq, fq, fastq.gz or fq.gz but file $fname doesn't have\n";}
	}
}



# tophat files:
# transcriptome index, gene annotation, bowtie index
my ($tophat_tr_index,$tophat_gtf,$tophat_bowtie_index);
if($genome eq "hg"){
	($tophat_tr_index,$tophat_gtf,$tophat_bowtie_index)=("/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/","/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/cuffcmp.combined.gtf","/users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/hg19/hg19");
}
if($genome eq "mm"){
	($tophat_tr_index,$tophat_gtf,$tophat_bowtie_index)=("/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/ENSEMBL_mm10_GRVm30/","/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/ENSEMBL_mm10_GRVm30/cuffcmp.combined.gtf","/users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/mm10/mm10");
}
if($genome eq "dr"){
	($tophat_tr_index,$tophat_gtf,$tophat_bowtie_index)=("/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/danRer10_ENSEMBL/","/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/danRer10_ENSEMBL/cuffcmp.combined.gtf","/users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/dr10/dr10");
}




#########
#########
my $call;
print "*************\nSANJUAN\n*************\n\n";
unless($test_run){
	print "Call:\nsanjuan -g $genome -c $conf -i $IRM -s $SuppJun -p $phred_code -l $library_type -a $adapter -g1 $g1_shortname -f1 ".join(" ",@g1_files)." -g2 $g2_shortname -f2 ".join(" ",@g2_files)." -o $output_dir -b $start_with -r $low_seq_req\n";
}else{
	print "Call:\nsanjuan -g $genome -c $conf -i $IRM -s $SuppJun -p $phred_code -l $library_type -a $adapter -g1 $g1_shortname -f1 ".join(" ",@g1_files)." -g2 $g2_shortname -f2 ".join(" ",@g2_files)." -o $output_dir -b $start_with -r $low_seq_req -t\n";
}
print "\n\n";

unless (-d $output_dir){print `mkdir -p $output_dir`;}
# here go all output and error messages
unless (-d "$output_dir/log_files"){print `mkdir -p $output_dir/log_files`;}

my $ret="Job ids:";
if($start_with ne "B"){
	# 1. triming and mapping
	# trim_galore needs python
	print "\n\n*************\nTrimming & Mapping\n*************\n\n";
	$call="perl $sanjuan_dir/preProcess_and_Map.pl $output_dir $genome $adapter $phred_code $library_type $start_with $test_run $tophat_tr_index $tophat_gtf $tophat_bowtie_index -g1 $g1_shortname @g1_files -g2 $g2_shortname @g2_files";
	$ret=`$call`;
	print $ret."\n";
}else{
	print "\n\n\nTrimming\n#####################\n\nSkipped\n\n\n\n\nMapping\n#####################\n\nSkipped\n\n\n\n\nMerging BAM files\n#####################\n\nSkipped\n\n";
}

my $job_ids=-1;
my @fs=split("\n",$ret);
if(@fs>0){
	chomp($fs[@fs-1]);
	$job_ids=$1 if $fs[@fs-1]=~/Job ids:(.+)/;
}

print "job ids from preProcess_and_Map: $job_ids\n\n\n";


my $merged_bam_file_1 = ($start_with eq "B")? $bam1 : $output_dir . '/TOPHAT_' . $g1_shortname."/accepted_hits_merged.bam";
my $merged_bam_file_2 = ($start_with eq "B")? $bam2 : $output_dir . '/TOPHAT_' . $g2_shortname."/accepted_hits_merged.bam";

print "merged_bam_file_1=$merged_bam_file_1\nmerged_bam_file_2=$merged_bam_file_2\n\n";

print "\n\n*************\nSplicing Analysis\n*************\n\n";
$call="perl $sanjuan_dir/SANJUAN_wrapper.pl $genome $RNAseq $conf $IRM $SuppJun $g1_shortname $g2_shortname $merged_bam_file_1 $merged_bam_file_2 $output_dir $job_ids $low_seq_req $test_run $sanjuan_dir $sanjuan_perllib $sanjuan_genomic_data_dir";
system($call);

exit(0);