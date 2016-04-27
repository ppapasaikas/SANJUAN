#!/usr/bin/perl
use strict;
use warnings;
my $Lqueue="long-sl65";		#Long queue
my $Squeue="short-sl65";	#Short queue
#Specify base directory for Output:
my $basedir=$ARGV[0];
my $genome=$ARGV[1]; 	#Specify species genome: hg -> human, mm-> mouse, dr-> zebrafish

# remove trailing "/" if exists
$basedir =~ s/\/$//;

my $adapter_seq=$ARGV[2];

#NRS	my $phred_code=$ARGV[3];

my $library_type=$ARGV[3];	#CinS 1S|1U|2S|2U

my $tpm=$ARGV[4];	#NinS. Two pass mapping

#NRS my $start_with = $ARGV[5];  #CinS M: mapping S: splicing 

my $test_run=$ARGV[5]; # 1 test run (statements are printed but not sent to cluster), 0 no test run

my $run_without_qsub=$ARGV[6]; # 1 = local run without qsub, 0 = with qsub

my $sanjuan_dir=$ARGV[7]; # location of SANJUAN program files

my $STAR_index=$ARGV[8];

#NRS my $tophat_bowtie_index=$ARGV[10];
my $N_processes=$ARGV[9];

my $phred_code=$ARGV[10];

#NRS	my $mate_inner_dist=$ARGV[12];
#NRS	my $mate_inner_dist_std_dev=$ARGV[13];
##########################################################################################################
##########################################################################################################
##########################################################################################################
# array of all ids of jobs which have been started within this script
my @all_job_ids=();
my ($job_id,$ret);
my $call;
# 0 -> no jobs would be started (test run) or have been started
# 1 -> at least one job would be started (test run) or have been started
my $jobs_started=0;
my $cond1_name="";
my $cond2_name="";
my @READ1=();
my @READ2=();
my @READS=(); #NinS single end sequencing files
my @READC=(); #NinS is set to either @READ1 or @READS depending on library type
my @COND=();
my $group="";
#NRS my $last_was_read2_file=1;
my @g1_files=();
my @g2_files=();
my @skipping_ok=();


for(my $i=11; $i<@ARGV; $i++){
	if(substr($ARGV[$i],0,3) eq "-g1"){
		$i++;
		$cond1_name=$ARGV[$i];
		$group=1;
		next;
	}

	if(substr($ARGV[$i],0,3) eq "-g2"){
		$i++;
		$cond2_name=$ARGV[$i];
		$group=2;
		next;
	}

	if($library_type=~/2/ && $ARGV[$i]=~/1\.(fq$|fastq$|fq\.gz$|fastq\.gz$|fq\.bz2$|fastq\.bz2$)/i){	#CinS
		push(@READ1,$ARGV[$i]);
		#NRS $last_was_read2_file=0;
		if($group==1){push(@g1_files,$ARGV[$i]);push(@COND,$cond1_name)};
		if($group==2){push(@g2_files,$ARGV[$i]);push(@COND,$cond2_name)};
		next;
	}

	if($library_type=~/2/ && $ARGV[$i]=~/2\.(fq$|fastq$|fq\.gz$|fastq\.gz$|fq\.bz2$|fastq\.bz2$)/i){	#CinS
		push(@READ2,$ARGV[$i]);
		#NRS $last_was_read2_file=1;
		if($group==1){push(@g1_files,$ARGV[$i])};
		if($group==2){push(@g2_files,$ARGV[$i])};
		next;
	}
	#NinS: pass single end sequencing files:
	if($library_type=~/1/ && $ARGV[$i]=~/\.(fq$|fastq$|fq\.gz$|fastq\.gz$|fq\.bz2$|fastq\.bz2$)/i){		#CinS
		push(@READS,$ARGV[$i]);
		if($group==1){push(@g1_files,$ARGV[$i]);push(@COND,$cond1_name)};
		if($group==2){push(@g2_files,$ARGV[$i]);push(@COND,$cond2_name)};
		next;
	}
}

@READC=($library_type=~/2/)? @READ1 : @READS;



for(my $i=0; $i<@READC; $skipping_ok[$i++]=1){}
my %skipping_ok;
$skipping_ok{$cond1_name}=1;
$skipping_ok{$cond2_name}=1;

print "\nCall of sub routine:\n";
print "preProcess_and_Map.pl $basedir $genome $adapter_seq $library_type $tpm $test_run $run_without_qsub $sanjuan_dir $STAR_index $N_processes -g1 $cond1_name ".join(",",@g1_files)." -g2 $cond2_name ".join(",",@g2_files)."\n";


########################################	M A P P I N G		##########################################

print `mkdir -p $basedir/MAPPING`;

print "\n\n\nMapping\n#####################\n\n";
my %bam_files=();
my %bam_files_count=();
my $readCommand="";
#my ($fq1,$fq2,$top_out);
my ($fq, $star_out);
my $trimCommand=($adapter_seq=~/\w+/)? "--clip3pAdapterSeq $adapter_seq" : "";

$,="\t";


my $outQSconversionAdd=0; # phred33
if($phred_code eq "phred64"){$outQSconversionAdd=-31;}

for my $cf (0..$#READC){
	my $cond=$COND[$cf];
	my $jname="MAP_" . $cond ."_$cf";
	#$fq1=$READ1[$cf];
	#$fq2=$READ2[$cf];
	$fq=($library_type=~/2/)? "$READ1[$cf] $READ2[$cf]" : $READS[$cf];
	
	if ($fq=~/\.gz$/){	#command for reading fastq files
	$readCommand="--readFilesCommand zcat";
	} elsif ($fq=~/\.bz2$/){
	$readCommand="--readFilesCommand bzip2 -c";	
	}	

	my $star_out=$basedir . '/MAPPING/' . $cond . "_" . "$cf" . "_";
	my $bampath=$star_out . "Aligned.sortedByCoord.out.bam";
	unless($bam_files{$cond}){$bam_files{$cond}=$bampath;
	} else {$bam_files{$cond}.=" $bampath";}
	$bam_files_count{$cond}++;

	# found already existing output 
	if(-f "$bampath" && $skipping_ok[$cf]){print "Skipping mapping of $fq because of already existing output $bampath\n\n"; next;}

	$skipping_ok{$cond}=0;
	$jobs_started=1;
	if($run_without_qsub==0){
		$call = "qsub -q $Lqueue -V -cwd -N $jname -o $basedir/log_files/02_out_mapping_${cf}.txt -e $basedir/log_files/02_err_mapping_${cf}.txt -pe smp $N_processes -l virtual_free=64G -l h_rt=48:00:00 -b y STAR --outQSconversionAdd $outQSconversionAdd --runThreadN $N_processes --genomeDir $STAR_index --readFilesIn $fq $readCommand $trimCommand --outSJfilterReads Unique --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 6 --alignSJDBoverhangMin 3 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --seedSearchStartLmax 50 --twopassMode $tpm --outFileNamePrefix $star_out";}
		else{
		$call = "STAR --outQSconversionAdd $outQSconversionAdd --runThreadN $N_processes --genomeDir $STAR_index --readFilesIn $fq $readCommand $trimCommand --outSJfilterReads Unique --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 6 --alignSJDBoverhangMin 3 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --seedSearchStartLmax 50 --twopassMode $tpm --outFileNamePrefix $star_out";
	}


	print "$call\n";
	if(!$test_run){
		$ret= `$call`;
		print "$ret\n\n";
		if($run_without_qsub==0){
   			($job_id)=$ret=~/job (\d+?) \(/;
			push(@all_job_ids,$job_id);
		}else{
			push(@all_job_ids,-3);
		}
	}

}



print "\n\n\nMerging BAM files\n#####################\n\n";

my $job_ids=join(",",@all_job_ids);
my $tmpdir="";
my $c=0;
foreach my $cond_name ($cond1_name,$cond2_name){
	my $merged_bam_file=$basedir . '/MAPPING/' . $cond_name . "_Aligned.sortedByCoord.out.merged.bam";
	$c++;

	if(-f $merged_bam_file && $skipping_ok{$cond_name}){print "Skipping merging of BAM files for condition $cond_name because of already existing output $merged_bam_file\n\n";next;}
	
	$jobs_started=1;

	# only one file
	if($bam_files_count{$cond_name}==1){
		if($run_without_qsub==0){
			$call="qsub -N SANJUAN_ln_s -q $Squeue -o $basedir/log_files/03_out_merging_bams_${c}.txt -e $basedir/log_files/03_err_merging_bams_${c}.txt -hold_jid $job_ids -V -cwd -l virtual_free=1G -l h_rt=00:20:00 $sanjuan_dir/ln_s_wrapper.sh $bam_files{$cond_name} $merged_bam_file";
		}else{
			$call="$sanjuan_dir/ln_s_wrapper.sh $bam_files{$cond_name} $merged_bam_file 1>$basedir/log_files/03_out_merging_bams_${c}.txt 2>$basedir/log_files/03_err_merging_bams_${c}.txt";
		}
	}else{
		if($run_without_qsub==0){
			$call = "qsub -N SANJUAN_BAMmerge_${cond_name}_${c} -q $Lqueue -o $basedir/log_files/03_out_merging_bams_${c}.txt -e $basedir/log_files/03_err_merging_bams_${c}.txt -hold_jid $job_ids -V -cwd -l virtual_free=40G -l h_rt=24:00:00 -b y samtools merge $merged_bam_file $bam_files{$cond_name}";
		}else{
			$call = "samtools merge $merged_bam_file $bam_files{$cond_name} 1>$basedir/log_files/03_out_merging_bams_${c}.txt 2>$basedir/log_files/03_err_merging_bams_${c}.txt";
		}
	}
	print "$call\n";
	if(!$test_run){
		$ret=`$call`;
		print "$ret\n\n";
		if($run_without_qsub==0){
			($job_id)=$ret=~/job (\d+?) \(/;
			push(@all_job_ids,$job_id);
		}else{
			push(@all_job_ids,-3);
		}
	}
}


# jobs would have been started but were not started because of test run
# we encode this by a job id -2 -> to be used in SANJUAN_wrapper.pl
if ($jobs_started && @all_job_ids==0){$all_job_ids[0]=-2;}
print "\n\n";
print "Job ids:".join(",",@all_job_ids);
