#!/usr/bin/perl
use strict;
use warnings;

### *S*plicing *AN*alysis & *JU*nction *AN*notation
### Software Requirements:
### samtools, bedtools, overlapSelect
### Perl Scripts: Found in SANJUAN directory
### Annotation Files: Found in SANJUAN directory/annotation files
### Ensembl_Transcript Junctions 'path2SANJUAN/annotation_files/...Transcript_Junctions.txt', EnsemblID2Name,
### Ensembl_Transcripts 'path2SANJUAN/annotation_files/...Transcripts.bed';

################### Subs ##############################################################################################################
#######################################################################################################################################
sub qsub { # $_[0]: command / job call as string; $_[1] is reference to all_job_ids_array; $_[2] output file of job; $_[3] reference to array @ENFORCE_RUN; $_[4] test_run=1 -> Yes 0-> real run 
		
	# ENFORCE_RUN or output does not exist
	if(${$_[3]}[0] || ! -f $_[2]){
		
		# delete output file (which might exist if run is enforced)
		if($_[4]){ # if test run
			print "deleting $_[2]\n";
		}else{
			unlink($_[2]);
		}
		
		print "$_[0]\n";
		
		# if not test run
		unless($_[4]){
			# delete output- and error-files
			unlink($1) if ($_[0] =~ /\-o (.+?) /);
			unlink($1) if ($_[0] =~ /\-e (.+?) /);
			sleep(1);
			my $ret = `$_[0]`;
			print "$ret\n\n";
			my ($job_id)=$ret=~/job (\d+?) \(/;
			push(@{$_[1]},$job_id);
		}
		
		# set $ENFORCE_RUN[0]=1 in order to enforce all later qsub calls to send the job to the cluster if their output file exists or not
		${$_[3]}[0]=1;
	
	}else{
		print "$_[0]\n... is skipped because of already existing output $_[2].\n\n";
	}
	print "\n";
}

# a reference to this array is a parameter to the sub qsub
# internally qsub checks if output of current job already exists; if yes, the job is not startet
# the output file is a parameter to qsub
# exception: when $ENFORCE_RUN[0]=1, 1. the output file is delete (if it exists), and 2. the job is sent to the cluster
# $ENFORCE_RUN[0] is set to 1 internally in qsub whenever it happens that one output file does not exists
# this enforces that all later runs of qsub are enforced to send the jobs to the cluster (if their output exists or not) 
my @ENFORCE_RUN=(0);

################# Setting Parameters ##################################################################################################
#######################################################################################################################################
my $genome=$ARGV[0];	#Specify species genome: hg-> human, mm-> mouse, dr-> zebrafish
my $RNAseq=$ARGV[1];	#Stranded "S" (e.g firstrand) or unstranded "U" RNAseq experiment. Default "S". 
my $conf=$ARGV[2];	#Specify stringency for Differentially Spliced Junctions (VeryHighConfidence -> VHC, HighConfidence -> HC, MediumConfidence -> MC)
my $IRM=$ARGV[3];	#IRM mode: Perform  High Sensitivity Intron Retention Anlaysis. Default 'N' 
my $SuppJun=$ARGV[4];	#Require Supporting Junction Evidence for IR identification (IRM mode). Default 'N'
 
#### Condition IDs ####
my $COND1=$ARGV[5];
my $COND2=$ARGV[6];

#### Location of (merged) bam files for COND1, COND2:
my $bam1=$ARGV[7];
my $bam2=$ARGV[8];

my $output_dir=$ARGV[9];
# remove trailing "/" if exists
$output_dir =~ s/\/$//;

# all previous jobs for whose ends we have to wait
my $job_ids=$ARGV[10];
my @all_job_ids=split(",",$job_ids);
my $job_id="";

# if job_ids="-1" -> no job has been started in preprocess_and_Map.pl
# if job_ids="175171,175172" -> jobs has been started in preprocess_and_Map.pl -> and therefor we have to run all jobs in SANJUAN_wrapper.pl
# if job_ids="-2" -> jobs would have been started in preprocess_and_Map.pl if it wasn't a test run -> and therefor we have to run all jobs in SANJUAN_wrapper.pl
if($all_job_ids[0] != -1){
	@ENFORCE_RUN=(1);
}

my $low_seq_req=$ARGV[11];
my $test_run=$ARGV[12]; # 1 test run (statements are printed but not sent to cluster), 0 no test run

# main paths variables
my $sanjuan_dir=$ARGV[13];
my $sanjuan_perllib=$ARGV[14]; 
my $sanjuan_genomic_data_dir=$ARGV[15];

my $prefix="SANJUAN";

print "Call of sub routine:\n";
print "SANJUAN_wrapper.pl $genome $RNAseq $conf $IRM $SuppJun $COND1 $COND2 $bam1 $bam2 $output_dir $job_ids $low_seq_req $test_run $sanjuan_dir $sanjuan_perllib $sanjuan_genomic_data_dir\n";


my ($ret,$genomePath,$Tx_bed,$ENS_Tx_Junc,$ENSid2Name);

#### Location of SANJUAN directory
if($genome eq "hg"){
	#### Location of bedtools genome files (originally from ~/ppapasaikas/SOFTWARE/bedtools2-2.20.1/genomes/ )
	$genomePath="$sanjuan_genomic_data_dir/genomes/human.hg19.genome";
	#### Location of Transcripts bed files (downloaded from UCSC table.browser | cut -f 1-6 OR gtf2Txbed.pl -used for zebrafish-)
	$Tx_bed="$sanjuan_genomic_data_dir/annotation_files/hg19_UCSC_Ensembl_Transcripts.bed";
	#### Location of Transcript Junctions files (downloaded gtf from UCSC table.browser | perl TxBed_J2Tx_from_gtf.pl)
	$ENS_Tx_Junc="$sanjuan_genomic_data_dir/annotation_files/hg19_UCSC_Ensembl_Transcript_Junctions.txt";
	#### Location of Transcript ID2 GeneName files (downloaded from UCSC table.browser)
	$ENSid2Name="$sanjuan_genomic_data_dir/annotation_files/hg19_UCSC_EnsemblTxID2name.txt";
}
if($genome eq "mm"){
	$genomePath="$sanjuan_genomic_data_dir/genomes/mouse.mm10.genome";
	$Tx_bed="$sanjuan_genomic_data_dir/annotation_files/mm10_UCSC_Ensembl_Transcripts.bed";
	$ENS_Tx_Junc="$sanjuan_genomic_data_dir/annotation_files/mm10_UCSC_Ensembl_Transcript_Junctions.txt";
	$ENSid2Name="$sanjuan_genomic_data_dir/annotation_files/mm10_UCSC_EnsemblTxID2name.txt";
}
if($genome eq "dr"){
	$genomePath="$sanjuan_genomic_data_dir/genomes/zebrafish.dr10.genome";
	$Tx_bed="$sanjuan_genomic_data_dir/annotation_files/dr10_EnsemblGRCz10_Transcripts.bed";
	$ENS_Tx_Junc="$sanjuan_genomic_data_dir/annotation_files/dr10_EnsemblGRCz10_Transcript_Junctions.txt";
	$ENSid2Name="$sanjuan_genomic_data_dir/annotation_files/dr10_EnsemblGRCz10TxID2Name.txt";
}

#### For debugging
#### Specify Skipping of preProcessing/Processing Steps ####
my %skip=();
$skip{1}=0;	#skip preprocessing step: -> bam file preprocessing
$skip{2}=0; #skip preprocessing step: -> Building of Non-Junction bam files
$skip{3}=0;	#skip: build SAM junction files, parsing, merging, slopping
$skip{4}=0;	#Skip  Run overlapSelect to find Neighboring Junctions and Junction-subsuming Transcripts
$skip{5}=0;	#skip: Calculate Differential Junction Efficiencies
$skip{6}=0;	#skip: Generate Intronic Segments and sorting
$skip{7}=0;	#skip: Calculate Coverage of Intronic Segments by NonJunction Reads:
$skip{8}=0;	#skip: Calculate Differential Intron Retention  
$skip{9}=0;	#skip: Annotate Differential Junctions
$skip{10}=0;#skip: Calculate Differential Intron Retention (IRM mode)
$skip{11}=0;#skip: Annotate Differential Introns

my $skip_nsteps=0;	#Number of processing steps to skip in analysis. e.g set to 2 to repeat analysis with different stringency
for (0..$skip_nsteps){ $skip{$_}=1; }

################### Main ##############################################################################################################
#######################################################################################################################################

#SANITIZE i.e remove unmapped reads and secondary alignments, keep only proper mates
my $PPoutbam1=$output_dir."/".$COND1 . '_ProperPReads.bam';
my $PPoutbam2=$output_dir."/".$COND2 . '_ProperPReads.bam';
my $Rmd_PPoutbam1=$output_dir."/".$COND1 . '_ProperPReads_rmdup.bam';
my $Rmd_PPoutbam2=$output_dir."/".$COND2 . '_ProperPReads_rmdup.bam';
my $NSort_PPoutbam1=$output_dir."/".$COND1 . '_ProperPReads_Nsorted';
my $NSort_PPoutbam2=$output_dir."/".$COND2 . '_ProperPReads_Nsorted';
my $paired = ($low_seq_req eq "N")? "-bf 0x2" : "-b";	#set to '-b' to keep non properly aligned mates / '-bf 0x2' to discard them
my $rmdup=0;	#set to 0 to keep PCR duplicates / 1 to discard them

unless ($skip{1}){
	print "\n\n\nRemoval of unmapped reads and secondary alignments\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_buildPPbam1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -e $output_dir/log_files/04_err_samtools_view_grp1.txt -o $PPoutbam1 -b y samtools view $paired -F 260 $bam1",\@all_job_ids,$PPoutbam1,\@ENFORCE_RUN,$test_run);	
	qsub("qsub -N ${prefix}_buildPPbam2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -e $output_dir/log_files/04_err_samtools_view_grp2.txt -o $PPoutbam2 -b y samtools view $paired -F 260 $bam2",\@all_job_ids,$PPoutbam2,\@ENFORCE_RUN,$test_run);

	if ($rmdup>0){
		$job_ids=join(",",@all_job_ids);
		qsub("qsub -N ${prefix}_rmdup1 -hold_jid $job_ids -V -cwd -pe -l virtual_free=32G -o $output_dir/log_files/04_out_samtools_rmdup_grp1.txt -e $output_dir/log_files/04_err_samtools_rmdup_grp1.txt -b y samtools rmdup $PPoutbam1 $Rmd_PPoutbam1",\@all_job_ids,$Rmd_PPoutbam1,\@ENFORCE_RUN,$test_run);
		qsub("qsub -N ${prefix}_rmdup2 -hold_jid $job_ids -V -cwd -pe -l virtual_free=32G -o $output_dir/log_files/04_out_samtools_rmdup_grp2.txt -e $output_dir/log_files/04_err_samtools_rmdup_grp2.txt -b y samtools rmdup $PPoutbam2 $Rmd_PPoutbam2",\@all_job_ids,$Rmd_PPoutbam2,\@ENFORCE_RUN,$test_run);
		
		$PPoutbam1=$Rmd_PPoutbam1;
		$PPoutbam2=$Rmd_PPoutbam2;
	}
}

### Build Non-Junction-Reads bam files.  (Keep Only pairs overlapping Junction Introns for both mates (pairtobed -type both))
### Generate bed and fix orientation of reads: Strand for mate /1 needs to be switched OR use -S in bedtools intersect
my $NJoutbam1=$output_dir."/".$COND1 . '_NoJunctReads.bam';
my $NJoutbam2=$output_dir."/".$COND2 . '_NoJunctReads.bam';
my $tempNJoutbed1=$output_dir."/".$COND1 . '_NoJunctReads_temp.bed';
my $tempNJoutbed2=$output_dir."/".$COND2 . '_NoJunctReads_temp.bed';
my $NJoutbed1=$output_dir."/".$COND1 . '_NoJunctReads.bed';
my $NJoutbed2=$output_dir."/".$COND2 . '_NoJunctReads.bed';

unless ($skip{2}){

	print "\n\n\nBuilding bed files for non-junction reads\n#####################\n\n";
	qsub("qsub -N $prefix -q short-sl65 -V -cwd -l virtual_free=1G -l h_rt=00:10:00 -o $output_dir/log_files/05_out_write_buildNJ1.txt -e $output_dir/log_files/05_err_write_buildNJ1.txt -b y perl $sanjuan_dir/job1.pl $PPoutbam1 $output_dir/buildNJ1.sh",\@all_job_ids,$NJoutbam1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N $prefix -q short-sl65 -V -cwd -l virtual_free=1G -l h_rt=00:10:00 -o $output_dir/log_files/05_out_write_buildNJ2.txt -e $output_dir/log_files/05_err_write_buildNJ2.txt -b y perl $sanjuan_dir/job1.pl $PPoutbam2 $output_dir/buildNJ2.sh",\@all_job_ids,$NJoutbam2,\@ENFORCE_RUN,$test_run);
	
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_buildNJbam1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $NJoutbam1 -e $output_dir/log_files/05_err_run_buildNJ1.txt -b y source $output_dir/buildNJ1.sh",\@all_job_ids,$NJoutbam1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_buildNJbam2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $NJoutbam2 -e $output_dir/log_files/05_err_run_buildNJ2.txt -b y source $output_dir/buildNJ2.sh",\@all_job_ids,$NJoutbam2,\@ENFORCE_RUN,$test_run);
	
	##### (Keep Only pairs overlapping Junction Introns for both mates (pairtobed -type both))
	### bedtools pairtobed -f 0.1 -type both -abam $NSort_PPoutbam1 -b JunctionIntrons.bed
	### bedtools pairtobed -f 0.1 -type both -abam $NSort_PPoutbam2 -b JunctionIntrons.bed

	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_buildNJbed1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $tempNJoutbed1 -e $output_dir/log_files/05_err_bamToBed_grp1.txt -b y bamToBed -i $NJoutbam1",\@all_job_ids,$tempNJoutbed1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_buildNJbed2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $tempNJoutbed2 -e $output_dir/log_files/05_err_bamToBed_grp2.txt -b y bamToBed -i $NJoutbam2",\@all_job_ids,$tempNJoutbed2,\@ENFORCE_RUN,$test_run);
	
	#####Fix orientation:
	$job_ids=join(",",@all_job_ids);
	my $tr='tr/\-\+/\+\-/ if \$_=~/\/1\s/';
	qsub("qsub -N ${prefix}_Fixbed1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -i $tempNJoutbed1 -o $NJoutbed1 -e $output_dir/log_files/05_err_fixOrient_grp1.txt -b y perl -pe \"\'$tr\'\"",\@all_job_ids,$NJoutbed1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_Fixbed1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -i $tempNJoutbed2 -o $NJoutbed2 -e $output_dir/log_files/05_err_fixOrient_grp2.txt -b y perl -pe \"\'$tr\'\"",\@all_job_ids,$NJoutbed2,\@ENFORCE_RUN,$test_run);
}


###################################################################################################################
my $JSAM1=$output_dir."/".$COND1 . '_JUNCT.sam';
my $JSAM2=$output_dir."/".$COND2 . '_JUNCT.sam';
my $Tophat_Junctions1=$output_dir."/"."Junctions_$COND1.bed";
my $Tophat_Junctions2=$output_dir."/"."Junctions_$COND2.bed";
my $SLOP1=$output_dir.'/Processed_pm' . '1000_Merged_Junctions.bed';
my $SLOP2=$output_dir.'/Processed_pm' . '10_Merged_Junctions.bed';

unless ($skip{3}){	#Get Junctions from SAM files, parse Cond1, Cond2, Merge and slop
	print "\n\n\nCreation of junction files\n#####################\n\n";
	#print `echo "samtools view $PPoutbam1 | awk $d6 ~ /N/' " > buildSAM1.sh`;
	#print `echo "samtools view $PPoutbam2 | awk $d6 ~ /N/' " > buildSAM2.sh`;
	
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N $prefix -q short-sl65 -V -cwd -l virtual_free=1G -l h_rt=00:10:00 -o $output_dir/log_files/06_out_write_buildSAM1.txt -e $output_dir/log_files/06_err_write_buildSAM1.txt -b y perl $sanjuan_dir/job2.pl $PPoutbam1 $output_dir/buildSAM1.sh",\@all_job_ids,$JSAM1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N $prefix -q short-sl65 -V -cwd -l virtual_free=1G -l h_rt=00:10:00 -o $output_dir/log_files/06_out_write_buildSAM2.txt -e $output_dir/log_files/06_err_write_buildSAM2.txt -b y perl $sanjuan_dir/job2.pl $PPoutbam2 $output_dir/buildSAM2.sh",\@all_job_ids,$JSAM2,\@ENFORCE_RUN,$test_run);
	
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_buildSAM1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $JSAM1 -e $output_dir/log_files/06_err_run_buildSAM1.txt -b y source $output_dir/buildSAM1.sh",\@all_job_ids,$JSAM1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_buildSAM2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $JSAM2 -e $output_dir/log_files/06_err_run_buildSAM1.txt -b y source $output_dir/buildSAM2.sh",\@all_job_ids,$JSAM2,\@ENFORCE_RUN,$test_run);
	
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_getJN1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $Tophat_Junctions1 -e $output_dir/log_files/06_err_get_juncts_grp1.txt -b y perl $sanjuan_dir/get_juncts.pl $JSAM1 $low_seq_req",\@all_job_ids,$Tophat_Junctions1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_getJN2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $Tophat_Junctions2 -e $output_dir/log_files/06_err_get_juncts_grp2.txt -b y perl $sanjuan_dir/get_juncts.pl $JSAM2 $low_seq_req",\@all_job_ids,$Tophat_Junctions2,\@ENFORCE_RUN,$test_run);

	$job_ids=join(",",@all_job_ids);
	qsub("qsub -q short-sl65 -N ${prefix}_merge_junctions -hold_jid $job_ids -V -cwd -l virtual_free=20G -o $output_dir/log_files/06_out_merge_juncts.txt -e $output_dir/log_files/06_err_merge_juncts.txt -b y perl $sanjuan_dir/merge_junctions.pl $Tophat_Junctions1 $Tophat_Junctions2 $output_dir/Merged_Junctions.bed",\@all_job_ids,"$output_dir/Merged_Junctions.bed",\@ENFORCE_RUN,$test_run);
	
	############### SLOP
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_SLOP1K -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $SLOP1 -e $output_dir/log_files/06_err_bedtools_slop1000.txt -b y bedtools slop -i $output_dir/Merged_Junctions.bed -g $genomePath -b 1000",\@all_job_ids,$SLOP1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_SLOP2K -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $SLOP2 -e $output_dir/log_files/06_err_bedtools_slop10.txt -b y bedtools slop -i $output_dir/Merged_Junctions.bed -g $genomePath -b 10",\@all_job_ids,$SLOP2,\@ENFORCE_RUN,$test_run);
}


#Run overlapSelect to find Neighboring Junctions and Junction-subsuming Transcripts
my $FileA=$SLOP1;
my $FileB=$SLOP2;
my $OUT_NJ_bed=$output_dir.'/olapSel_JUNCpm1000_JUNCpm10.bed';
my $FileA2=$Tx_bed;
my $FileB2=$output_dir.'/Merged_Junctions.bed';
my $OUT_Jun2Tx_bed=$output_dir.'/olapSel_ENSTX_JUNC.bed';
unless ($skip{4}){
	print "\n\n\nFinding neighboring junctions and junction-subsuming transcripts\n#####################\n\n";

	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_OLSLN -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/07_out_overlapSelect_1.txt -e $output_dir/log_files/07_err_overlapSelect_1.txt -b y overlapSelect -strand -mergeOutput $FileA $FileB $OUT_NJ_bed",\@all_job_ids,$OUT_NJ_bed,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_OLSLT -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/07_out_overlapSelect_2.txt -e $output_dir/log_files/07_err_overlapSelect_2.txt -b y overlapSelect -strand -overlapThreshold=1.0 -mergeOutput $FileA2 $FileB2 $OUT_Jun2Tx_bed",\@all_job_ids,$OUT_Jun2Tx_bed,\@ENFORCE_RUN,$test_run);
}


#Calculate Differential Junction Efficiencies:
my $Proc_Junctions1=$Tophat_Junctions1;
my $Proc_Junctions2=$Tophat_Junctions2;
my $olapSel_NJunc12=$OUT_NJ_bed;
my $olapSel_Junc2Tx=$OUT_Jun2Tx_bed;
my $OUT_calc_HC_JEFF=$output_dir."/Diff_Junctions_".$conf.".txt";
my $OUT_calc_LC_JEFF=$output_dir.'/Diff_Junctions_LC.txt';
my $OUT_calc_NC_JEFF=$output_dir.'/Diff_Junctions_NC.txt';
my @par=($Proc_Junctions1,$Proc_Junctions2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc);
unless ($skip{5}){
	print "\n\n\nCalculation of high- and low-confidence differential junctions\n#####################\n\n";
	
	$job_ids=join(",",@all_job_ids);
	# XXX should be ommit this file because now we have the output for all junctions with parameter NC?				
	qsub("qsub -N ${prefix}_LCDJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $OUT_calc_LC_JEFF -e $output_dir/log_files/08_err_calcJUNCTeff_Lconf.txt -b y perl $sanjuan_dir/calc_JUNCT_efficiency.pl $par[0] $par[1] $par[2] $par[3] LC $par[4] $sanjuan_perllib",\@all_job_ids,$OUT_calc_LC_JEFF,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_HCDJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $OUT_calc_HC_JEFF -e $output_dir/log_files/08_err_calcJUNCTeff_conf.txt -b y perl $sanjuan_dir/calc_JUNCT_efficiency.pl $par[0] $par[1] $par[2] $par[3] $conf $par[4] $sanjuan_perllib",\@all_job_ids,$OUT_calc_HC_JEFF,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_NCDJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $OUT_calc_NC_JEFF -e $output_dir/log_files/08_err_calcJUNCTeff_NC.txt -b y perl $sanjuan_dir/calc_JUNCT_efficiency.pl $par[0] $par[1] $par[2] $par[3] NC $par[4] $sanjuan_perllib",\@all_job_ids,$OUT_calc_NC_JEFF,\@ENFORCE_RUN,$test_run);
}


#Generate Intronic Segments
my $OUT_INTR_SEGM=$output_dir.'/Junctions_IntronicSegments.bed';
my $OUT_INTR_SEGM_SORTED=$output_dir.'/Junctions_IntronicSegments_sorted.bed';
unless ($skip{6}){
	print "\n\n\nGeneration of intronic segments bed file\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N $prefix -hold_jid $job_ids -V -cwd -l virtual_free=20G -o $OUT_INTR_SEGM -e $output_dir/log_files/09_err_juncts2Introns.txt -b y perl $sanjuan_dir/tophat_junctions2IntronSegments.pl $output_dir/Merged_Junctions.bed",\@all_job_ids,$OUT_INTR_SEGM,\@ENFORCE_RUN,$test_run);
		
	print "\n\n\nSorting intronic segments file\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	#print `sort -k1,1 -k2,2n $OUT_INTR_SEGM > $OUT_INTR_SEGM_SORTED`;
	qsub("qsub -N $prefix -hold_jid $job_ids -V -cwd -o $OUT_INTR_SEGM_SORTED -e $output_dir/log_files/09_err_sortIntronSegments.txt -l virtual_free=5G $sanjuan_dir/sort_wrapper.sh $OUT_INTR_SEGM",\@all_job_ids,$OUT_INTR_SEGM_SORTED,\@ENFORCE_RUN,$test_run);
}


#Calculate Coverage of Intronic Segments by NonJunction Reads:
my $NoJunctReads_bed1=$NJoutbed1;		#Generated during Preprocessing
my $NoJunctReads_bed2=$NJoutbed2;		#Generated during Preprocessing
my $OUT_Coverage_IntrSegm1= "$output_dir/$COND1"."_IntrSegm_coverage.bed";
my $OUT_Coverage_IntrSegm2= "$output_dir/$COND2"."_IntrSegm_coverage.bed";
my $str_spec="-s";
$str_spec="" if $RNAseq eq "U";
unless ($skip{7}){
	print "\n\n\nCalculation intronic segments read coverage\n#####################\n\n";
	
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_IntrCov1 -hold_jid $job_ids -V -cwd -pe ompi 2 -l virtual_free=32G -o $OUT_Coverage_IntrSegm1 -e $output_dir/log_files/10_err_bedtools_intersect_grp1.txt -b y bedtools intersect -c $str_spec -sorted -a $OUT_INTR_SEGM_SORTED -b $NoJunctReads_bed1",\@all_job_ids,$OUT_Coverage_IntrSegm1,\@ENFORCE_RUN,$test_run);
	qsub("qsub -N ${prefix}_IntrCov2 -hold_jid $job_ids -V -cwd -pe ompi 2 -l virtual_free=32G -o $OUT_Coverage_IntrSegm2 -e $output_dir/log_files/10_err_bedtools_intersect_grp1.txt -b y bedtools intersect -c $str_spec -sorted -a $OUT_INTR_SEGM_SORTED -b $NoJunctReads_bed2",\@all_job_ids,$OUT_Coverage_IntrSegm2,\@ENFORCE_RUN,$test_run);
}


#Calculate Differential Intron Retention
@par=($OUT_calc_HC_JEFF,$OUT_Coverage_IntrSegm1,$OUT_Coverage_IntrSegm2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc,$Proc_Junctions1,$Proc_Junctions2,"noIRM");
my $OUT_INTR_RET=$output_dir."/Diff_RetIntr.txt";
unless ($skip{8}){
	print "\n\n\nCalculation of differential intron retention\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_DIR -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $OUT_INTR_RET -e $output_dir/log_files/11_err_calcIntronRet.txt -b y perl $sanjuan_dir/calc_INTRON_retention.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] $par[7] $par[8] $sanjuan_perllib",\@all_job_ids,$OUT_INTR_RET,\@ENFORCE_RUN,$test_run);
}


#Annotate Differential Junctions
@par=($OUT_calc_HC_JEFF,$OUT_calc_LC_JEFF,$olapSel_Junc2Tx,$OUT_INTR_RET,$ENSid2Name,$ENS_Tx_Junc);
my $OUT_ANNOT=$output_dir."/Annotated_Diff_Junctions.txt";
unless ($skip{9}){
	print "\n\n\nAnnotation of differential junctions\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	qsub("qsub -N ${prefix}_ADJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $OUT_ANNOT -e $output_dir/log_files/12_err_annotateDiffJuncts.txt -b y perl $sanjuan_dir/annotate_Diff_Used_Junctions.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $conf $COND1 $COND2",\@all_job_ids,$OUT_ANNOT,\@ENFORCE_RUN,$test_run);
}


if ($IRM eq 'Y'){

	#Calculate Differential Intron Retention (IRM mode)
	@par=($OUT_calc_LC_JEFF,$OUT_Coverage_IntrSegm1,$OUT_Coverage_IntrSegm2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc,$Proc_Junctions1,$Proc_Junctions2,"IRM");
	my $OUT_IRM=$output_dir."/Diff_RetIntr_IRM.txt";
	unless ($skip{10}){
		print "\n\n\nCalculation of differential intron retention (IRM mode)\n#####################\n\n";
		$job_ids=join(",",@all_job_ids);
		qsub("qsub -N ${prefix}_IRM -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $OUT_IRM -e $output_dir/log_files/13_err_calcIntronRet_IRM.txt -b y perl $sanjuan_dir/calc_INTRON_retention.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] $par[7] $par[8] $sanjuan_perllib",\@all_job_ids,$OUT_IRM,\@ENFORCE_RUN,$test_run);
	}

	# should this be inside the if($IRM eq "Y") or not?
	#Annotate Differential Introns
	@par=($OUT_calc_LC_JEFF,$OUT_calc_LC_JEFF,$olapSel_Junc2Tx,$OUT_IRM,$ENSid2Name,$ENS_Tx_Junc,$SuppJun);
	my $OUT_IANNOT="Annotated_Diff_Introns.txt";
	unless ($skip{11}){
		print "\n\n\nAnnotation of differential introns\n#####################\n\n";
		$job_ids=join(",",@all_job_ids);
		qsub("qsub -N ${prefix}_ADI -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $OUT_IANNOT -e $output_dir/log_files/14_err_annotateDiffIntrons_IRM.txt -b y perl $sanjuan_dir/annotate_Diff_Used_Introns.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] VHC $COND1 $COND2",\@all_job_ids,$OUT_IANNOT,\@ENFORCE_RUN,$test_run);
	}
}

print "\n\n\n******************\nDone: ";
if($test_run){print "did not send jobs to cluster because this was a test run\n\n";}else{print "sent all jobs to cluster\n\n";}
print "all job ids: ".join(",",@all_job_ids)."\n";