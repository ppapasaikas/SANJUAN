#!/usr/bin/perl
use strict;
use warnings;

my $genome=$ARGV[1]; 	#Specify species genome: hg -> human, mm-> mouse, dr-> zebrafish

#Specify base directory for Output:
my $basedir=$ARGV[0];
# remove trailing "/" if exists
$basedir =~ s/\/$//;

my $adapter_seq=$ARGV[2];

my $phred_code=$ARGV[3];

my $library_type=$ARGV[4];

my $start_with = $ARGV[5];  # T: trimming, M: mapping

my $test_run=$ARGV[6]; # 1 test run (statements are printed but not sent to cluster), 0 no test run

my $run_without_qsub=$ARGV[7]; # 1 = local run without qsub, 0 = with qsub

my $sanjuan_dir=$ARGV[8]; # location of SANJUAN program files

my $tophat_tr_index=$ARGV[9];
my $tophat_bowtie_index=$ARGV[10];
my $N_processes=$ARGV[11];
my $mate_inner_dist=$ARGV[12];
my $mate_inner_dist_std_dev=$ARGV[13];
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
my @COND=();
my $group="";
my $last_was_read2_file=1;
my @g1_files=();
my @g2_files=();
my @skipping_ok=();

for(my $i=14; $i<@ARGV; $i++){
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

	if($last_was_read2_file){
		push(@READ1,$ARGV[$i]);
		$last_was_read2_file=0;
		if($group==1){push(@g1_files,$ARGV[$i]);push(@COND,$cond1_name)};
		if($group==2){push(@g2_files,$ARGV[$i]);push(@COND,$cond2_name)};
		next;
	}

	if(!$last_was_read2_file){
		push(@READ2,$ARGV[$i]);
		$last_was_read2_file=1;
		if($group==1){push(@g1_files,$ARGV[$i])};
		if($group==2){push(@g2_files,$ARGV[$i])};
		next;
	}
}

for(my $i=0; $i<@READ1; $skipping_ok[$i++]=1){}
my %skipping_ok;
$skipping_ok{$cond1_name}=1;
$skipping_ok{$cond2_name}=1;

print "Call of sub routine:\n";
print "preProcess_and_Map.pl $basedir $genome $adapter_seq $phred_code $library_type $start_with $test_run $run_without_qsub $sanjuan_dir $tophat_tr_index $tophat_bowtie_index $N_processes -g1 $cond1_name ".join(",",@g1_files)." -g2 $cond2_name ".join(",",@g2_files)."\n";

########################################	T R I M M I N G		##########################################
if($start_with eq "T"){
	#Specify adapter sequence to be trimmed. To find adapter: minion search-adapter -i FASTQFILE.gz
	my $ad=$adapter_seq;	#(Default adapter sequence for CRG facility, RNAseq )
	my $bc="";
	my $ad_bc=$ad . $bc;
	my $odir="";
	my $jname;
	
	print "\n\n\nTrimming\n#####################\n\n";
	# naming of output files? always fq.gz or depending on input also fastq.gz -> yes, trim galore names output files always fq.gz
	for my $cf (0..$#READ1){
		my $cond=$COND[$cf];
		$odir=$basedir . "/TRIM_$cond";
		print `mkdir -p $odir`;
		my $outfile1=$odir."/".$1."_val_1.fq.gz" if ("/".$READ1[$cf])=~ /\/(\w+)(\.fastq\.gz$|\.fastq\.gzip$|\.fq\.gz$|\.fq\.gzip$|\.fq$|\.fastq$)/;
		my $outfile2=$odir."/".$1."_val_2.fq.gz" if ("/".$READ2[$cf])=~ /\/(\w+)(\.fastq\.gz$|\.fastq\.gzip$|\.fq\.gz$|\.fq\.gzip$|\.fq$|\.fastq$)/;
		# if output is already there
		if(-f $outfile1 && -f $outfile2){
			print "Skipping trimming of files $READ1[$cf] and $READ2[$cf] because of already existing output $outfile1 and $outfile2.\n\n";
			next;
		}
		$skipping_ok[$cf]=0;
		$jname="TRIM_$cond"."_$cf";
		if($run_without_qsub==0){
			$call= "qsub -q long-sl65 -V -cwd -N $jname -o $basedir/log_files/01_out_trim_${cf}.txt -e $basedir/log_files/01_err_trim_${cf}.txt -l virtual_free=40G -l h_rt=48:00:00 $sanjuan_dir/trim_galore_wrapper.sh --$phred_code --gzip --stringency 3 -q 0 -a $ad_bc -a2 $ad_bc --length 19 --paired $READ1[$cf] $READ2[$cf] -o $odir";
		}else{
			$call= "$sanjuan_dir/trim_galore_wrapper.sh --$phred_code --gzip --stringency 3 -q 0 -a $ad_bc -a2 $ad_bc --length 19 --paired $READ1[$cf] $READ2[$cf] -o $odir 1>$basedir/log_files/01_out_trim_${cf}.txt 2>$basedir/log_files/01_err_trim_${cf}.txt";
		}
		print "$call\n";
		$jobs_started=1;
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
}else{
	print "\n\n\nTrimming\n#####################\n\nSkipped\n\n.";
	
}

########################################	M A P P I N G		##########################################
##TO GET BOWTIE INDEXES:
#print `wget -P /Volumes/HD2/RESEARCH/DATA/BOWTIE2_INDEXES/mm10/ ftp://ftp.cbcb.umd.edu/pub/data/bowtie2_indexes/incl/mm10.zip`;	#hg19.zip/mm10.zip
#OR from ilumina's iGENOMES
#print `wget -P /Volumes/HD2/RESEARCH/DATA/BOWTIE2_INDEXES/mm10/ ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz`;
#then: tar -zxvf ....gz;

print `mkdir -p $basedir/TOPHAT_$cond1_name`;
print `mkdir -p $basedir/TOPHAT_$cond2_name`;

print "\n\n\nMapping\n#####################\n\n";
my $job_ids=join(",",@all_job_ids);
my %bam_files=();
my %bam_files_count=();
my ($fq1,$fq2,$top_out);
for my $cf (0..$#READ1){
	my $cond=$COND[$cf];
	my $odir=$basedir . "/TRIM_$cond";
	my $jname="TOP_$cond"."_$cf";
	if($start_with eq "M"){# use input fastq files because they are already trimmed
		$fq1=$READ1[$cf];
		$fq2=$READ2[$cf];
	}else{ # use output from trimming
		$fq1=$odir."/".$1."_val_1.fq.gz" if ("/".$READ1[$cf])=~ /\/(\w+)(\.fastq\.gz$|\.fastq\.gzip$|\.fq\.gz$|\.fq\.gzip$|\.fq$|\.fastq$)/;
		$fq2=$odir."/".$1."_val_2.fq.gz" if ("/".$READ2[$cf])=~ /\/(\w+)(\.fastq\.gz$|\.fastq\.gzip$|\.fq\.gz$|\.fq\.gzip$|\.fq$|\.fastq$)/;
	}
	$top_out="$basedir/TOPHAT_$cond/$cf";

	unless($bam_files{$cond}){$bam_files{$cond}="$top_out/accepted_hits.bam";
	}else{$bam_files{$cond}.=" $top_out/accepted_hits.bam";}
	$bam_files_count{$cond}++;

	# found already existing output 
	if(-f "$top_out/accepted_hits.bam" && $skipping_ok[$cf]){print "Skipping mapping of $fq1 and $fq2 because of already existing output $top_out/accepted_hits.bam\n\n"; next;}

	$skipping_ok{$cond}=0;
	$jobs_started=1;
	if($run_without_qsub==0){
		$call = "qsub -q long-sl65 -V -cwd -N $jname -o $basedir/log_files/02_out_mapping_${cf}.txt -e $basedir/log_files/02_err_mapping_${cf}.txt -hold_jid $job_ids -pe smp $N_processes -l virtual_free=64G -l h_rt=48:00:00 -b y tophat2 --no-mixed --library-type $library_type -o $top_out --transcriptome-index=$tophat_tr_index -p $N_processes --inner-mate-dist $mate_inner_dist --mate-std-dev $mate_inner_dist_std_dev -i 50 -I 800000 -x 1 $tophat_bowtie_index $fq1 $fq2";
	}else{
		$call = "tophat2 --no-mixed --library-type $library_type -o $top_out --transcriptome-index=$tophat_tr_index -p $N_processes --mate-inner-dist $mate_inner_dist --mate-std-dev $mate_inner_dist_std_dev -i 50 -I 800000 -x 1 $tophat_bowtie_index $fq1 $fq2 1>$basedir/log_files/02_out_mapping_${cf}.txt 2>$basedir/log_files/02_err_mapping_${cf}.txt";
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

	#Only for CNAG sequencing
	##print `qsub -q long-sl65 -V -cwd -N $cond._TOP -hold_jid $job_ids -e ~/temp -pe ompi 12 -l virtual_free=64G -l h_rt=48:00:00 -b y tophat2 --no-mixed --library-type fr-unstranded -o $top_out --transcriptome-index=/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/ -G /users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/iGENOMES_UCSC_hg19_clean/cuffcmp.combined.gtf -p 12 -r 20 --mate-std-dev 45 -i 50 -I 800000 -x 1 /users/jvalcarcel/ppapasaikas/BOWTIE2_INDEXES/hg19/hg19 $fq1 $fq2`;
}


print "\n\n\nMerging BAM files\n#####################\n\n";

$job_ids=join(",",@all_job_ids);
my $tmpdir="";
my $c=0;
foreach my $cond_name ($cond1_name,$cond2_name){
	$tmpdir=$basedir . '/TOPHAT_' . $cond_name;
	my $merged_bam_file=$tmpdir."/accepted_hits_merged.bam";
	$c++;

	if(-f $merged_bam_file && $skipping_ok{$cond_name}){print "Skipping merging of BAM files for condition $cond_name because of already existing output $merged_bam_file\n\n";next;}
	
	$jobs_started=1;

	# only one file
	if($bam_files_count{$cond_name}==1){
		if($run_without_qsub==0){
			$call="qsub -N SANJUAN_ln_s -q short-sl65 -o $basedir/log_files/03_out_merging_bams_${c}.txt -e $basedir/log_files/03_err_merging_bams_${c}.txt -hold_jid $job_ids -V -cwd -l virtual_free=1G -l h_rt=00:20:00 ~/crg/projects/2015_sanjuan/dev/ln_s_wrapper.sh $bam_files{$cond_name} $merged_bam_file";
		}else{
			$call="~/crg/projects/2015_sanjuan/dev/ln_s_wrapper.sh $bam_files{$cond_name} $merged_bam_file 1>$basedir/log_files/03_out_merging_bams_${c}.txt 2>$basedir/log_files/03_err_merging_bams_${c}.txt";
		}
	}else{
		if($run_without_qsub==0){
			$call = "qsub -N SANJUAN_SAMmerge_${cond_name}_${c} -q long-sl65 -o $basedir/log_files/03_out_merging_bams_${c}.txt -e $basedir/log_files/03_err_merging_bams_${c}.txt -hold_jid $job_ids -V -cwd -l virtual_free=40G -l h_rt=24:00:00 -b y samtools merge $merged_bam_file $bam_files{$cond_name}";
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
