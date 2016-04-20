#!/usr/bin/perl
use strict;
use warnings;

### *S*plicing *AN*alysis & *JU*nction *AN*notation
### Software Requirements:
### samtools, bedtools, overlapSelect
### Perl Scripts: Found in SANJUAN lib directory
### Annotation Files: Found in SANJUAN db/SANJUAN_annotation_files:
### Transcripts.bed Transcript_Junctions.txt, TxID2Name,


################### Subs ##############################################################################################################
#######################################################################################################################################
sub run_cmd { # $_[0]: command / job call as string / qsub part; $_[1]: command / job / program call; $_[2] is reference to all_job_ids_array; $_[3] output file of job; $_[4] reference to array @ENFORCE_RUN; $_[5] test_run=1 -> Yes 0-> real run; $_[6]: run_without_qsub (0/1) 
		
	my $call="";
	# if use qsub
	if($_[6]==0){
		$call="$_[0] $_[1]";
	}else{
		$call=$_[1];
	}
		
	# ENFORCE_RUN or output does not exist
	if(${$_[4]}[0] || ! -f $_[3]){
		
		# delete output file (which might exist if run is enforced)
		if($_[5]){ # if test run
			print "deleting $_[3]\n";
		}else{
			unlink($_[3]);
		}

		print "$call\n";	

		# if not test run
		unless($_[5]){
			# delete output- and error-files
			unlink($1) if ($_[0] =~ /\-o (.+?) /);
			unlink($1) if ($_[0] =~ /\-e (.+?) /);
			sleep(1);
			my $ret = `$call`;
			print "$ret\n\n";
			# run with qsub
			if($_[6]==0){
				my ($job_id)=$ret=~/job (\d+?) \(/;
				push(@{$_[2]},$job_id);
			}
		}

		# set $ENFORCE_RUN[0]=1 in order to enforce all later qsub calls to send the job to the cluster if their output file exists or not
		${$_[4]}[0]=1;
	
	}else{
		print "$call\n... is skipped because of already existing output $_[3].\n\n";
	}
	print "\n";
}

# a reference to this array is a parameter to the sub qsub
# internally qsub checks if output of current job already exists; if yes, the job is not started
# the output file is a parameter to qsub
# exception: when $ENFORCE_RUN[0]=1, 1. the output file is deletes (if it exists), and 2. the job is sent to the cluster
# $ENFORCE_RUN[0] is set to 1 internally in qsub whenever it happens that one output file does not exists
# this enforces that all later runs of qsub are enforced to send the jobs to the cluster (if their output exists or not) 
my @ENFORCE_RUN=(0);

################# Setting Parameters ##################################################################################################
#######################################################################################################################################
my $genome=$ARGV[0];	#Specify species genome: hg-> human, mm-> mouse, dr-> zebrafish
my $library_type=$ARGV[1];	#1S|2S|1U|2U: Single end 1 or Pair-end 2 Stranded "S" (e.g firstrand) or unstranded "U" RNAseq experiment. Default "S". 
my $RNAseq=($library_type=~/S/)? "S":"U";

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
my $run_without_qsub=$ARGV[13]; # 0=> use qsub, 1=> don't use qsub
my $N_processes=$ARGV[14];
my $rmdup_as_argument=$ARGV[15];

# main paths variables
my $sanjuan_dir=$ARGV[16];
my $sanjuan_perllib=$ARGV[17]; 
my $sanjuan_genomic_data_dir=$ARGV[18];

my $prefix="SANJUAN";

print "Call of sub routine:\n";
print "SANJUAN_wrapper.pl $genome $RNAseq $conf $IRM $SuppJun $COND1 $COND2 $bam1 $bam2 $output_dir $job_ids $low_seq_req $test_run $run_without_qsub $N_processes $sanjuan_dir $sanjuan_perllib $sanjuan_genomic_data_dir\n";

my ($ret,$genomePath,$Tx_bed,$ENS_Tx_Junc,$ENSid2Name);
#### Location of bedtools genome files (originally from ~/ppapasaikas/SOFTWARE/bedtools2-2.20.1/genomes/ )
$genomePath="$sanjuan_genomic_data_dir/genomes/${genome}.genome";
### Location of Transcripts bed files (downloaded from UCSC table.browser | cut -f 1-6 OR gtf2Txbed.pl -used for zebrafish-)
$Tx_bed="$sanjuan_genomic_data_dir/SANJUAN_annotation_files/${genome}_Transcripts.bed";
### Location of Transcript Junctions files (downloaded gtf from UCSC table.browser | perl TxBed_J2Tx_from_gtf.pl)
$ENS_Tx_Junc="$sanjuan_genomic_data_dir/SANJUAN_annotation_files/${genome}_Transcript_Junctions.txt";
### Location of Transcript ID2 GeneName files (downloaded from UCSC table.browser)
$ENSid2Name="$sanjuan_genomic_data_dir/SANJUAN_annotation_files/${genome}_TxID2Name.txt";

#### For debugging
#### Specify Skipping of preProcessing/Processing Steps ####
my %skip=();
$skip{1}=0;	#skip preprocessing step: -> bam file preprocessing
#$skip{2}=0; 	#skip preprocessing step: -> Building of Non-Junction bam files. Deprecated. Now coverage of intronic segments is calculated drectly from BAM file
$skip{2}=0;	#skip: build SAM junction files, parsing, merging, slopping
$skip{3}=0;	#Skip  Run overlapSelect to find Neighboring Junctions and Junction-subsuming Transcripts
$skip{4}=0;	#skip: Calculate Differential Junction Efficiencies
$skip{5}=0;	#skip: Generate Intronic Segments and sorting
$skip{6}=0;	#skip: Calculate Coverage of Intronic Segments by NonJunction Reads:
$skip{7}=0;	#skip: Calculate Differential Intron Retention  
$skip{8}=0;	#skip: Annotate Differential Junctions
$skip{9}=0;#skip: Calculate Differential Intron Retention (IRM mode)
$skip{10}=0;#skip: Annotate Differential Introns

my $skip_nsteps=0;	#Number of processing steps to skip in analysis. e.g set to 2 to repeat analysis with different stringency
for (0..$skip_nsteps){ $skip{$_}=1; }

################### Main ##############################################################################################################
#######################################################################################################################################

#SANITIZE i.e remove unmapped reads and secondary alignments, keep only proper mates. Deprecated after STAR transition.
#Now only removal of duplicate reads...
#NRS my $PPoutbam1=$output_dir."/".$COND1 . '_ProperPReads.bam';
#NRS my $PPoutbam2=$output_dir."/".$COND2 . '_ProperPReads.bam';
my $Rmd_bam1=$output_dir."/".$COND1 . '_rmdup.bam';
my $Rmd_bam2=$output_dir."/".$COND2 . '_rmdup.bam';
#NRS my $paired = ($low_seq_req eq "N")? "-bf 0x2" : "-b";	#set to '-b' to keep non properly aligned mates / '-bf 0x2' to discard them
my $rmdup=$rmdup_as_argument;	#set to 0 to keep PCR duplicates / 1 to discard them

unless ($skip{1} || $rmdup==0){
	#print "\n\n\nRemoval of unmapped reads and secondary alignments\n#####################\n\n";
	print "\n\n\nRemoval of duplicate reads\n#####################\n\n";
#NRS	$job_ids=join(",",@all_job_ids);
#NRS	run_cmd("qsub -N ${prefix}_buildPPbam1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -e $output_dir/log_files/04_err_samtools_view_grp1.txt -o $output_dir/log_files/04_out_samtools_view_grp1.txt -b y","samtools view $paired -F 260 -o $PPoutbam1 $bam1",\@all_job_ids,$PPoutbam1,\@ENFORCE_RUN,$test_run,$run_without_qsub);	
#NRS	run_cmd("qsub -N ${prefix}_buildPPbam2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -e $output_dir/log_files/04_err_samtools_view_grp2.txt -o $output_dir/log_files/04_out_samtools_view_grp2.txt -b y","samtools view $paired -F 260 -o $PPoutbam2 $bam2",\@all_job_ids,$PPoutbam2,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	if ($library_type=~/1/){
		$job_ids=join(",",@all_job_ids);
		run_cmd("qsub -N ${prefix}_rmdup1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/05_out_samtools_rmdup_grp1.txt -e $output_dir/log_files/05_err_samtools_rmdup_grp1.txt -b y","samtools rmdup -s $bam1 $Rmd_bam1",\@all_job_ids,$Rmd_bam1,\@ENFORCE_RUN,$test_run,$run_without_qsub);
		run_cmd("qsub -N ${prefix}_rmdup2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/05_out_samtools_rmdup_grp2.txt -e $output_dir/log_files/05_err_samtools_rmdup_grp2.txt -b y","samtools rmdup -s $bam2 $Rmd_bam2",\@all_job_ids,$Rmd_bam2,\@ENFORCE_RUN,$test_run,$run_without_qsub);

	}
	else 	{
		$job_ids=join(",",@all_job_ids);
		run_cmd("qsub -N ${prefix}_rmdup1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/05_out_samtools_rmdup_grp1.txt -e $output_dir/log_files/05_err_samtools_rmdup_grp1.txt -b y","samtools rmdup $bam1 $Rmd_bam1",\@all_job_ids,$Rmd_bam1,\@ENFORCE_RUN,$test_run,$run_without_qsub);
		run_cmd("qsub -N ${prefix}_rmdup2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/05_out_samtools_rmdup_grp2.txt -e $output_dir/log_files/05_err_samtools_rmdup_grp2.txt -b y","samtools rmdup $bam2 $Rmd_bam2",\@all_job_ids,$Rmd_bam2,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	}
$bam1=$Rmd_bam1;
$bam2=$Rmd_bam2;
}



###################################################################################################################
my $JSAM1=$output_dir."/".$COND1 . '_JUNCT.sam';
my $JSAM2=$output_dir."/".$COND2 . '_JUNCT.sam';
my $Mapped_Junctions1=$output_dir."/"."Junctions_$COND1.bed";
my $Mapped_Junctions2=$output_dir."/"."Junctions_$COND2.bed";
my $SLOP1=$output_dir.'/Processed_pm' . '1000_Merged_Junctions.bed';
my $SLOP2=$output_dir.'/Processed_pm' . '10_Merged_Junctions.bed';

unless ($skip{2}){	#Get Junctions from SAM files, parse Cond1, Cond2, Merge and slop
	print "\n\n\nCreation of junction files\n#####################\n\n";
	#print `echo "samtools view $PPoutbam1 | awk $d6 ~ /N/' " > buildSAM1.sh`;
	#print `echo "samtools view $PPoutbam2 | awk $d6 ~ /N/' " > buildSAM2.sh`;
	
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N $prefix -q short-sl65 -V -cwd -l virtual_free=1G -l h_rt=00:10:00 -o $output_dir/log_files/08_out_write_buildSAM1.txt -e $output_dir/log_files/08_err_write_buildSAM1.txt -b y","perl $sanjuan_dir/job2.pl $bam1 $output_dir/buildSAM1.sh $JSAM1",\@all_job_ids,"$output_dir/buildSAM1.sh",\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N $prefix -q short-sl65 -V -cwd -l virtual_free=1G -l h_rt=00:10:00 -o $output_dir/log_files/08_out_write_buildSAM2.txt -e $output_dir/log_files/08_err_write_buildSAM2.txt -b y","perl $sanjuan_dir/job2.pl $bam2 $output_dir/buildSAM2.sh $JSAM2",\@all_job_ids,"$output_dir/buildSAM2.sh",\@ENFORCE_RUN,$test_run,$run_without_qsub);
	
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N ${prefix}_buildSAM1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/09_out_run_buildSAM1.txt -e $output_dir/log_files/09_err_run_buildSAM1.txt -b y","$output_dir/buildSAM1.sh",\@all_job_ids,$JSAM1,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N ${prefix}_buildSAM2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/09_out_run_buildSAM2.txt -e $output_dir/log_files/09_err_run_buildSAM2.txt -b y","$output_dir/buildSAM2.sh",\@all_job_ids,$JSAM2,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N ${prefix}_getJN1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/10_out_get_juncts_grp1.txt -e $output_dir/log_files/10_err_get_juncts_grp1.txt -b y","perl $sanjuan_dir/get_juncts.pl $JSAM1 $low_seq_req $Mapped_Junctions1",\@all_job_ids,$Mapped_Junctions1,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N ${prefix}_getJN2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/10_out_get_juncts_grp2.txt -e $output_dir/log_files/10_err_get_juncts_grp2.txt -b y","perl $sanjuan_dir/get_juncts.pl $JSAM2 $low_seq_req $Mapped_Junctions2",\@all_job_ids,$Mapped_Junctions2,\@ENFORCE_RUN,$test_run,$run_without_qsub);

	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -q short-sl65 -N ${prefix}_merge_junctions -hold_jid $job_ids -V -cwd -l virtual_free=20G -o $output_dir/log_files/11_out_merge_juncts.txt -e $output_dir/log_files/11_err_merge_juncts.txt -b y","perl $sanjuan_dir/merge_junctions.pl $Mapped_Junctions1 $Mapped_Junctions2 $output_dir/Merged_Junctions.bed",\@all_job_ids,"$output_dir/Merged_Junctions.bed",\@ENFORCE_RUN,$test_run,$run_without_qsub);
	
	############### SLOP
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N ${prefix}_SLOP1K -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/12_out_bedtools_slop1000.txt -e $output_dir/log_files/12_err_bedtools_slop1000.txt -b y","$sanjuan_dir/bedtools_slop.sh $output_dir/Merged_Junctions.bed $genomePath 1000 $SLOP1",\@all_job_ids,$SLOP1,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N ${prefix}_SLOP2K -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/12_out_bedtools_slop10.txt -e $output_dir/log_files/12_err_bedtools_slop10.txt -b y","$sanjuan_dir/bedtools_slop.sh $output_dir/Merged_Junctions.bed $genomePath 10 $SLOP2",\@all_job_ids,$SLOP2,\@ENFORCE_RUN,$test_run,$run_without_qsub);
}


#Run overlapSelect to find Neighboring Junctions and Junction-subsuming Transcripts
my $FileA=$SLOP1;
my $FileB=$SLOP2;
my $OUT_NJ_bed=$output_dir.'/olapSel_JUNCpm1000_JUNCpm10.bed';
my $FileA2=$Tx_bed;
my $FileB2=$output_dir.'/Merged_Junctions.bed';
my $OUT_Jun2Tx_bed=$output_dir.'/olapSel_ENSTX_JUNC.bed';
unless ($skip{3}){
	print "\n\n\nFinding neighboring junctions and junction-subsuming transcripts\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N ${prefix}_INTSN -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/13_out_bedtools_intersect_NJ.txt -e $output_dir/log_files/13_err_bedtools_intersect_NJ.txt -b y","$sanjuan_dir/bedtools_intersect_NJ.sh $FileB $FileA $OUT_NJ_bed",\@all_job_ids,$OUT_NJ_bed,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N ${prefix}_INTST -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/13_out_bedtools_intersect_ST.txt -e $output_dir/log_files/13_err_bedtools_intersect_ST.txt -b y","$sanjuan_dir/bedtools_intersect_ST.sh $FileB2 $FileA2 $OUT_Jun2Tx_bed",\@all_job_ids,$OUT_Jun2Tx_bed,\@ENFORCE_RUN,$test_run,$run_without_qsub);
}


#Calculate Differential Junction Efficiencies:
my $Proc_Junctions1=$Mapped_Junctions1;
my $Proc_Junctions2=$Mapped_Junctions2;
my $olapSel_NJunc12=$OUT_NJ_bed;
my $olapSel_Junc2Tx=$OUT_Jun2Tx_bed;
my $OUT_calc_HC_JEFF=$output_dir."/Diff_Junctions_".$conf.".txt";
my $OUT_calc_LC_JEFF=$output_dir.'/Diff_Junctions_LC.txt';
my $OUT_calc_NC_JEFF=$output_dir.'/Diff_Junctions_NC.txt';
my @par=($Proc_Junctions1,$Proc_Junctions2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc);
unless ($skip{4}){
	print "\n\n\nCalculation of high- and low-confidence differential junctions\n#####################\n\n";
	
	$job_ids=join(",",@all_job_ids);
	# XXX should this file be ommit because now we have the output for all junctions with parameter NC?				
	run_cmd("qsub -N ${prefix}_LCDJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/14_out_calcJUNCTeff_LC.txt -e $output_dir/log_files/14_err_calcJUNCTeff_LC.txt -b y","perl $sanjuan_dir/calc_JUNCT_efficiency.pl $par[0] $par[1] $par[2] $par[3] LC $par[4] $sanjuan_perllib $OUT_calc_LC_JEFF",\@all_job_ids,$OUT_calc_LC_JEFF,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N ${prefix}_HCDJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/14_out_calcJUNCTeff_HC.txt -e $output_dir/log_files/14_err_calcJUNCTeff_HC.txt -b y","perl $sanjuan_dir/calc_JUNCT_efficiency.pl $par[0] $par[1] $par[2] $par[3] $conf $par[4] $sanjuan_perllib $OUT_calc_HC_JEFF",\@all_job_ids,$OUT_calc_HC_JEFF,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N ${prefix}_NCDJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/14_out_calcJUNCTeff_NC.txt -e $output_dir/log_files/14_err_calcJUNCTeff_NC.txt -b y","perl $sanjuan_dir/calc_JUNCT_efficiency.pl $par[0] $par[1] $par[2] $par[3] NC $par[4] $sanjuan_perllib $OUT_calc_NC_JEFF",\@all_job_ids,$OUT_calc_NC_JEFF,\@ENFORCE_RUN,$test_run,$run_without_qsub);
}


#Generate Intronic Segments
my $OUT_INTR_SEGM=$output_dir.'/Junctions_IntronicSegments.bed';
my $OUT_INTR_SEGM_SORTED=$output_dir.'/Junctions_IntronicSegments_sorted.bed';
unless ($skip{5}){
	print "\n\n\nGeneration of intronic segments bed file\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N $prefix -hold_jid $job_ids -V -cwd -l virtual_free=20G -o $output_dir/log_files/15_out_juncts2Introns.txt -e $output_dir/log_files/15_err_juncts2Introns.txt -b y","perl $sanjuan_dir/Mapped_junctions2IntronSegments.pl $output_dir/Merged_Junctions.bed $OUT_INTR_SEGM",\@all_job_ids,$OUT_INTR_SEGM,\@ENFORCE_RUN,$test_run,$run_without_qsub);
		
	print "\n\n\nSorting intronic segments file\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	#print `sort -k1,1 -k2,2n $OUT_INTR_SEGM > $OUT_INTR_SEGM_SORTED`;
	run_cmd("qsub -N $prefix -hold_jid $job_ids -V -cwd -o $output_dir/log_files/16_out_sortIntronSegments.txt -e $output_dir/log_files/16_err_sortIntronSegments.txt -l virtual_free=5G","$sanjuan_dir/sort_wrapper.sh $OUT_INTR_SEGM $OUT_INTR_SEGM_SORTED",\@all_job_ids,$OUT_INTR_SEGM_SORTED,\@ENFORCE_RUN,$test_run,$run_without_qsub);
}


#Calculate Coverage of Intronic Segments by NonJunction Reads:
my $OUT_Coverage_IntrSegm1= "$output_dir/$COND1"."_IntrSegm_coverage.bed";
my $OUT_Coverage_IntrSegm2= "$output_dir/$COND2"."_IntrSegm_coverage.bed";
unless ($skip{6}){
	print "\n\n\nCalculation of intronic segments read coverage\n#####################\n\n";
	
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N ${prefix}_IntrCov1 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/17_out_bedtools_intersect_grp1.txt -e $output_dir/log_files/17_err_bedtools_intersect_grp1.txt -b y","$sanjuan_dir/bedtools_intersect.sh $RNAseq $OUT_INTR_SEGM_SORTED $bam1 $OUT_Coverage_IntrSegm1",\@all_job_ids,$OUT_Coverage_IntrSegm1,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	run_cmd("qsub -N ${prefix}_IntrCov2 -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/17_out_bedtools_intersect_grp2.txt -e $output_dir/log_files/17_err_bedtools_intersect_grp2.txt -b y","$sanjuan_dir/bedtools_intersect.sh $RNAseq $OUT_INTR_SEGM_SORTED $bam2 $OUT_Coverage_IntrSegm2",\@all_job_ids,$OUT_Coverage_IntrSegm2,\@ENFORCE_RUN,$test_run,$run_without_qsub);
}


#Calculate Differential Intron Retention
@par=($OUT_calc_HC_JEFF,$OUT_Coverage_IntrSegm1,$OUT_Coverage_IntrSegm2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc,$Proc_Junctions1,$Proc_Junctions2,"noIRM");
my $OUT_INTR_RET=$output_dir."/Diff_RetIntr.txt";
unless ($skip{7}){
	print "\n\n\nCalculation of differential intron retention\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N ${prefix}_DIR -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/18_out_calcIntronRet.txt -e $output_dir/log_files/18_err_calcIntronRet.txt -b y","perl $sanjuan_dir/calc_INTRON_retention.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] $par[7] $par[8] $sanjuan_perllib $OUT_INTR_RET",\@all_job_ids,$OUT_INTR_RET,\@ENFORCE_RUN,$test_run,$run_without_qsub);
}

#Annotate Differential Junctions
@par=($OUT_calc_HC_JEFF,$OUT_calc_LC_JEFF,$olapSel_Junc2Tx,$OUT_INTR_RET,$ENSid2Name,$ENS_Tx_Junc);
my $OUT_ANNOT=$output_dir."/Annotated_Diff_Junctions.txt";
unless ($skip{8}){
	print "\n\n\nAnnotation of differential junctions\n#####################\n\n";
	$job_ids=join(",",@all_job_ids);
	run_cmd("qsub -N ${prefix}_ADJ -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/19_out_annotateDiffJuncts.txt -e $output_dir/log_files/19_err_annotateDiffJuncts.txt -b y","perl $sanjuan_dir/annotate_Diff_Used_Junctions.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $conf $COND1 $COND2 $OUT_ANNOT",\@all_job_ids,$OUT_ANNOT,\@ENFORCE_RUN,$test_run,$run_without_qsub);
}

if ($IRM eq 'Y'){
	#Calculate Differential Intron Retention (IRM mode)
	@par=($OUT_calc_LC_JEFF,$OUT_Coverage_IntrSegm1,$OUT_Coverage_IntrSegm2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc,$Proc_Junctions1,$Proc_Junctions2,"IRM");
	my $OUT_IRM=$output_dir."/Diff_RetIntr_IRM.txt";
	unless ($skip{10}){
		print "\n\n\nCalculation of differential intron retention (IRM mode)\n#####################\n\n";
		$job_ids=join(",",@all_job_ids);
		run_cmd("qsub -N ${prefix}_IRM -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/20_out_calcIntronRet_IRM.txt -e $output_dir/log_files/20_err_calcIntronRet_IRM.txt -b y","perl $sanjuan_dir/calc_INTRON_retention.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] $par[7] $par[8] $sanjuan_perllib $OUT_IRM",\@all_job_ids,$OUT_IRM,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	}

	# should this be inside the if($IRM eq "Y") or not?
	#Annotate Differential Introns
	@par=($OUT_calc_LC_JEFF,$OUT_calc_LC_JEFF,$olapSel_Junc2Tx,$OUT_IRM,$ENSid2Name,$ENS_Tx_Junc,$SuppJun);
	my $OUT_IANNOT=$output_dir."/Annotated_Diff_Introns.txt";
	unless ($skip{11}){
		print "\n\n\nAnnotation of differential introns\n#####################\n\n";
		$job_ids=join(",",@all_job_ids);
		run_cmd("qsub -N ${prefix}_ADI -hold_jid $job_ids -V -cwd -l virtual_free=32G -o $output_dir/log_files/21_out_annotateDiffIntrons_IRM.txt -e $output_dir/log_files/21_err_annotateDiffIntrons_IRM.txt -b y","perl $sanjuan_dir/annotate_Diff_Used_Introns.pl $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] VHC $COND1 $COND2 $OUT_IANNOT",\@all_job_ids,$OUT_IANNOT,\@ENFORCE_RUN,$test_run,$run_without_qsub);
	}
}

print "\n\n\n******************\nDone: ";
if($test_run){
	if($run_without_qsub==0){
		print "did not send jobs to cluster because this was a test run\n\n";
	}else{
		print "did not start jobs because this was a test run\n\n";
	}
}else{
	if($run_without_qsub==0){
		print "sent all jobs to cluster\n\n";
	}else{
		print "ran all jobs\n\n";
	}
}

print "all job ids: ".join(",",@all_job_ids)."\n";
