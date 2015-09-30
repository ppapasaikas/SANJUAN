### *S*plicing *AN*alysis & *JU*nction *AN*notation
### Software Requirements:
### samtools, bedtools, overlapSelect
### Perl Scripts: Found in SANJUAN directory
### Annotation Files: Found in SANJUAN directory/annotation files
### Ensembl_Transcript Junctions 'path2SANJUAN/annotation_files/...Transcript_Junctions.txt', EnsemblID2Name,
### Ensembl_Transcripts 'path2SANJUAN/annotation_files/...Transcripts.bed';

$genome="hg";	#Specify species genome: hg-> human, mm-> mouse, dr-> zebrafish
$RNAseq="S";	#Stranded "S" (e.g firstrand) or unstranded "U" RNAseq experiment. Default "S". 
$conf='HC';	#Specify stringency for Differentially Spliced Junctions (VeryHighConfidence -> VHC, HighConfidence -> HC, MediumConfidence -> MC)
$IRM='Y';	#IRM mode: Perform  High Sensitivity Intron Retention Anlaysis. Default 'N' 
$SuppJun='N';	#Require Supporting Junction Evidence for IR identification (IRM mode). Default 'N'
 
#### Condition IDs ####
$COND1="CNT";
$COND2="SSA";

$prefix="LuVi"; # Cluster Jobs IDs Prefix, this should be in the form: LETTERS_

#### Location of (merged) bam files for COND1, COND2:
$bam1='/users/jvalcarcel/ppapasaikas/LUISA/STAR_Control_2ndM/Aligned.sortedByCoord.out.bam';
$bam2='/users/jvalcarcel/ppapasaikas/LUISA/STAR_SSA_2ndM/Aligned.sortedByCoord.out.bam';



#######################################################################################################################################
#######################################################################################################################################

#### Location of SANJUAN directory
$perlPath='/users/jvalcarcel/ppapasaikas/SANJUAN_DEV/';
#### Location of bedtools genome files (originally from ~/ppapasaikas/SOFTWARE/bedtools2-2.20.1/genomes/ )
$genomePath=$perlPath . 'genomes/human.hg19.genome' if $genome eq 'hg';
$genomePath=$perlPath . 'genomes/mouse.mm10.genome' if $genome eq 'mm';
$genomePath=$perlPath . 'genomes/zebrafish.dr10.genome' if $genome eq 'dr';
#### Location of Transcripts bed files (downloaded from UCSC table.browser | cut -f 1-6 OR gtf2Txbed.pl -used for zebrafish-)
$Tx_bed=$perlPath . 'annotation_files/hg19_UCSC_Ensembl_Transcripts.bed' if $genome eq 'hg';
$Tx_bed=$perlPath . 'annotation_files/mm10_UCSC_Ensembl_Transcripts.bed' if $genome eq 'mm';
$Tx_bed=$perlPath . 'annotation_files/dr10_EnsemblGRCz10_Transcripts.bed' if $genome eq 'dr';
#### Location of Transcript Junctions files (downloaded gtf from UCSC table.browser | perl TxBed_J2Tx_from_gtf.pl)
$ENS_Tx_Junc=$perlPath . 'annotation_files/hg19_UCSC_Ensembl_Transcript_Junctions.txt' if $genome eq 'hg';
$ENS_Tx_Junc=$perlPath . 'annotation_files/mm10_UCSC_Ensembl_Transcript_Junctions.txt' if $genome eq 'mm';
$ENS_Tx_Junc=$perlPath . 'annotation_files/dr10_EnsemblGRCz10_Transcript_Junctions.txt' if $genome eq 'dr';
#### Location of Transcript ID2 GeneName files (downloaded from UCSC table.browser)
$ENSid2Name=$perlPath . 'annotation_files/hg19_UCSC_EnsemblTxID2name.txt' if $genome eq 'hg';
$ENSid2Name=$perlPath . 'annotation_files/mm10_UCSC_EnsemblTxID2name.txt' if $genome eq 'mm';
$ENSid2Name=$perlPath . 'annotation_files/dr10_EnsemblGRCz10TxID2Name.txt' if $genome eq 'dr';


#### Specify Skipping of preProcessing/Processing Steps ####
$preskip{1}=0;	#Skip preprocessing step: -> bam file preprocessing
$preskip{2}=0; 	#Skip preprocessing step: -> Building of Non-Junction bam files

$skip{1}=0;	#Skip individual step -> build SAM junction files, parsing, merging, slopping
$skip{2}=0;	#Skip individual step -> overlapSelect for Juntions/Transcripts
$skip{3}=0;	#Skip individual step -> Differential Junction Efficiency Calculation
$skip{4}=0;	#Skip individual step -> Intronic Segment Generation and sorting
$skip{5}=0;	#Skip individual step -> Intronic Segment Coverage
$skip{6}=0;	#Skip individual step -> Differential Intron Retention Calculation
$skip{7}=0;	#Skip individual step -> Annotation File

$skip_nsteps=0;	#Number of processing steps to skip in analysis. e.g set to 2 to repeat analysis with different stringency

for (0..$skip_nsteps){
$skip{$_}=1;
}
$preskip{1}=1 if $skip_nsteps>0;
$preskip{2}=1 if $skip_nsteps>0;

print "\n";

$Lconf='LC';
$Lconf='HC' if $conf eq 'IRM';

########################### TIME CONSUMING PREPROCESSING (Cluster Jobs) ##################################

#SANITIZE i.e remove unmapped reads and secondary alignments, keep only proper mates
$PPoutbam1=$COND1 . '_ProperPReads.bam';
$PPoutbam2=$COND2 . '_ProperPReads.bam';
$Rmd_PPoutbam1=$COND1 . '_ProperPReads_rmdup.bam';
$Rmd_PPoutbam2=$COND2 . '_ProperPReads_rmdup.bam';
$NSort_PPoutbam1=$COND1 . '_ProperPReads_Nsorted';
$NSort_PPoutbam2=$COND2 . '_ProperPReads_Nsorted';
$paired='-bf 0x2';	#set to '-b' to keep non properly aligned mates / '-bf 0x2' to discard them
$rmdup=0;		#set to 0 to keep PCR duplicates / 1 to discard them
unless ($preskip{1}){
print `rm  -f *_ProperPReads*`;
print `qsub -N $prefix.buildPPbam1 -V -cwd -pe ompi 1 -l virtual_free=32G -o $PPoutbam1 -b y samtools view $paired -F 260 $bam1 >> flush.tmp`;
print `qsub -N $prefix.buildPPbam2 -V -cwd -pe ompi 1 -l virtual_free=32G -o $PPoutbam2 -b y samtools view $paired -F 260 $bam2 >> flush.tmp`;
	if ($rmdup>0){
	print `qsub -N $prefix.rmdup1 -hold_jid $prefix.buildPPbam1 -V -cwd -pe ompi 1 -l virtual_free=32G -b y samtools rmdup $PPoutbam1 $Rmd_PPoutbam1 >> flush.tmp`;
	print `qsub -N $prefix.rmdup2 -hold_jid $prefix.buildPPbam2 -V -cwd -pe ompi 1 -l virtual_free=32G -b y samtools rmdup $PPoutbam2 $Rmd_PPoutbam2 >> flush.tmp`;
	$PPoutbam1=$Rmd_PPoutbam1;
	$PPoutbam2=$Rmd_PPoutbam2;
	}

#print `qsub -N $prefix.SortPPbam1 -hold_jid $prefix.buildPPbam1,$prefix.rmdup1 -V -cwd -pe ompi 1 -l virtual_free=32G -b y samtools sort -n $PPoutbam1 $NSort_PPoutbam1 >> flush.tmp`;
#print `qsub -N $prefix.SortPPbam2 -hold_jid $prefix.buildPPbam2,$prefix.rdmup2 -V -cwd -pe ompi 1 -l virtual_free=32G -b y samtools sort -n $PPoutbam2 $NSort_PPoutbam2 >> flush.tmp`;
print `qsub -N $prefix.SortDummy -cwd -V -hold_jid $prefix.SortPPbam1,$prefix.SortPPbam2,$prefix.buildPPbam1,$prefix.rmdup1,$prefix.buildPPbam2,$prefix.rmdup2 -sync y -b y echo >> flush.tmp`;
print `rm -f flush.tmp`;
print `rm -f $prefix.*`;
}

if ($rmdup>0){
$PPoutbam1=$Rmd_PPoutbam1;
$PPoutbam2=$Rmd_PPoutbam2;
}



### Build Non-Junction-Reads bam files.  (Keep Only pairs overlapping Junction Introns for both mates (pairtobed -type both))
### Generate bed and fix orientation of reads: Strand for mate /1 needs to be switched OR use -S in bedtools intersect
$NJoutbam1=$COND1 . '_NoJunctReads.bam';
$NJoutbam2=$COND2 . '_NoJunctReads.bam';
$tempNJoutbed1=$COND1 . '_NoJunctReads_temp.bed';
$tempNJoutbed2=$COND2 . '_NoJunctReads_temp.bed';
$NJoutbed1=$COND1 . '_NoJunctReads.bed';
$NJoutbed2=$COND2 . '_NoJunctReads.bed';
unless ($preskip{2}){
$d6='\'\$6';
$d1='\$1';
$tr='tr/\-\+/\+\-/ if \$_=~/\/1\s/';
print `rm  -f *_NoJunctReads*`;
print `echo "samtools view -h $PPoutbam1 | awk $d6 !~ /N/ || $d1 ~ /^@/' | samtools view -bS - " > buildNJ1.sh`;
print `echo "samtools view -h $PPoutbam2 | awk $d6 !~ /N/ || $d1 ~ /^@/' | samtools view -bS - " > buildNJ2.sh`;
print `qsub -N $prefix.buildNJbam1 -V -cwd -pe ompi 1 -l virtual_free=32G -o $NJoutbam1 -b y source buildNJ1.sh > flush.tmp`;
print `qsub -N $prefix.buildNJbam2 -V -cwd -pe ompi 1 -l virtual_free=32G -o $NJoutbam2 -b y source buildNJ2.sh >> flush.tmp`;

##### (Keep Only pairs overlapping Junction Introns for both mates (pairtobed -type both))
### bedtools pairtobed -f 0.1 -type both -abam $NSort_PPoutbam1 -b JunctionIntrons.bed
### bedtools pairtobed -f 0.1 -type both -abam $NSort_PPoutbam2 -b JunctionIntrons.bed

print "Building bed files for Non-Junction Reads...\n";
print `qsub -N $prefix.buildNJbed1 -hold_jid $prefix.buildNJbam1 -V -cwd -pe ompi 1 -l virtual_free=32G -o $tempNJoutbed1 -b y bamToBed -i $NJoutbam1 >> flush.tmp`;
print `qsub -N $prefix.buildNJbed2 -hold_jid $prefix.buildNJbam2 -V -cwd -pe ompi 1 -l virtual_free=32G -o $tempNJoutbed2 -b y bamToBed -i $NJoutbam2 >> flush.tmp`;
#####Fix orientation:
print `qsub -N $prefix.Fixbed1 -hold_jid $prefix.buildNJbed1 -V -cwd -pe ompi 1 -l virtual_free=32G -i $tempNJoutbed1 -o $NJoutbed1 -b y perl -pe \"\'$tr\'\" >> flush.tmp`;
print `qsub -N $prefix.Fixbed2 -hold_jid $prefix.buildNJbed2 -V -cwd -pe ompi 1 -l virtual_free=32G -i $tempNJoutbed2 -o $NJoutbed2 -b y perl -pe \"\'$tr\'\" >> flush.tmp`;
print `rm -f $prefix.*`;
}



unless ($preskip{2}){
print `qsub -N $prefix.FixDummy -cwd -V -hold_jid $prefix.Fixbed1,$prefix.Fixbed2 -sync y -b y echo >> flush.tmp`;
print `rm -f flush.tmp`;
print `rm -f *_temp.bed`;
print "...Done building bed files for Non-Junction Reads\n";
print `rm -f $prefix.*`;
}







###################################################################################################################
$JSAM1=$COND1 . '_JUNCT.sam';
$JSAM2=$COND2 . '_JUNCT.sam';
$Tophat_Junctions1="Junctions_$COND1.bed";
$Tophat_Junctions2="Junctions_$COND2.bed";
$SLOP1='Processed_pm' . '1000_Merged_Junctions.bed';
$SLOP2='Processed_pm' . '10_Merged_Junctions.bed';

unless ($skip{1}){	#Get Junctions from SAM files, parse Cond1, Cond2, Merge and slop 
print "...Creating Junction Files\n";
$d6='\'\$6';
print `rm -f $JSAM1 $JSAM2 $Tophat_Junctions1 $Tophat_Junctions2`;
print `echo "samtools view $PPoutbam1 | awk $d6 ~ /N/' " > buildSAM1.sh`;
print `echo "samtools view $PPoutbam2 | awk $d6 ~ /N/' " > buildSAM2.sh`;
print `qsub -N $prefix.buildSAM1 -V -cwd -pe ompi 1 -l virtual_free=32G -o $JSAM1 -b y source buildSAM1.sh > flush.tmp`;
print `qsub -N $prefix.buildSAM2 -V -cwd -pe ompi 1 -l virtual_free=32G -o $JSAM2 -b y source buildSAM2.sh >> flush.tmp`;
print `qsub -N $prefix.SAMDummy -cwd -V -hold_jid $prefix.buildSAM1,$prefix.buildSAM2 -sync y -b y echo >> flush.tmp`;

print `qsub -N $prefix.getJN1 -V -cwd -pe ompi 1 -l virtual_free=32G -o $Tophat_Junctions1 -b y perl $perlPath\'get_juncts.pl' $JSAM1 > flush.tmp`;
print `qsub -N $prefix.getJN2 -V -cwd -pe ompi 1 -l virtual_free=32G -o $Tophat_Junctions2 -b y perl $perlPath\'get_juncts.pl' $JSAM2 >> flush.tmp`;
print `qsub -N $prefix.getJNDummy -cwd -V -hold_jid $prefix.getJN1,$prefix.getJN2 -sync y -b y echo >> flush.tmp`;

print `rm $prefix.*`;


################ Merge Junction Files
open (IN1, $Tophat_Junctions1);
		while (<IN1>){
		$line=$_;
		chomp $line;
		@mat=split /\t/,$line;
		$count{$mat[3]}=$mat[4];
		}
	close IN1;
	open (IN2, $Tophat_Junctions2);
		while (<IN2>){
		$line=$_;
		chomp $line;
		@mat=split /\t/,$line;
		$count{$mat[3]}+=$mat[4];
		}
	close IN2;
	open (OUT, ">Merged_Junctions.bed");
		foreach $id (keys %count){
		@INF=split /_/,$id;
		print OUT "$INF[0]\t$INF[1]\t$INF[2]\t$id\t$count{$id}\t$INF[3]\n";#$fcount{$id}\t$flowcells\n";
		}
	close OUT;

############### SLOP
print `rm -f $SLOP1 $SLOP2`;
print `qsub -N $prefix.SLOP1K -V -cwd -pe ompi 1 -l virtual_free=32G -o $SLOP1 -b y bedtools slop -i Merged_Junctions.bed -g $genomePath -b 1000 > flush.tmp`;
print `qsub -N $prefix.SLOP10 -V -cwd -pe ompi 1 -l virtual_free=32G -o $SLOP2 -b y bedtools slop -i Merged_Junctions.bed -g $genomePath -b 10 >> flush.tmp`;
print `qsub -N $prefix.SAMDummy -cwd -V -hold_jid $prefix.SLOP1K,$prefix.SLOP10 -sync y -b y echo >> flush.tmp`;
print `rm $prefix.*`;
}





#Run overlapSelect to find Neighboring Junctions and Junction-subsuming Transcripts
$FileA=$SLOP1;
$FileB=$SLOP2;
$OUT_NJ_bed='olapSel_JUNCpm1000_JUNCpm10.bed';
$FileA2=$Tx_bed;
$FileB2='Merged_Junctions.bed';
$OUT_Jun2Tx_bed='olapSel_ENSTX_JUNC.bed';
unless ($skip{2}){
print `rm -f $OUT_NJ_bed $OUT_Jun2Tx_bed`;
print "...Finding Neighboring Junctions and junction-subsuming Transcripts\n";
print `qsub -N $prefix.OLSLN -V -cwd -pe ompi 1 -l virtual_free=32G -b y overlapSelect -strand -mergeOutput $FileA $FileB $OUT_NJ_bed  > flush.tmp`;
print `qsub -N $prefix.OLSLT -V -cwd -pe ompi 1 -l virtual_free=32G -b y overlapSelect -strand -overlapThreshold=1.0 -mergeOutput $FileA2 $FileB2 $OUT_Jun2Tx_bed  >> flush.tmp`;
print `qsub -N $prefix.OLSLDummy -cwd -V -hold_jid $prefix.OLSLN,$prefix.OLSLT -sync y -b y echo >> flush.tmp`;
print `rm -f $prefix*`;
print "...Done\n";
}




#Calculate Differential Junction Efficiencies:
$Proc_Junctions1=$Tophat_Junctions1;
$Proc_Junctions2=$Tophat_Junctions2;
$olapSel_NJunc12=$OUT_NJ_bed;
$olapSel_Junc2Tx=$OUT_Jun2Tx_bed;
$OUT_calc_HC_JEFF='Diff_Junctions_HC.txt';
$OUT_calc_LC_JEFF='Diff_Junctions_LC.txt';
@par=($Proc_Junctions1,$Proc_Junctions2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc);
unless ($skip{3}){
print `rm -f $OUT_calc_HC_JEFF $OUT_calc_LC_JEFF`;
print "...Calculating High and Low Confidence Differential Junctions\n";
print `qsub -N $prefix.LCDJ -V -cwd -pe ompi 1 -l virtual_free=32G -o $OUT_calc_LC_JEFF -b y perl $perlPath\'calc_JUNCT_efficiency.pl' $par[0] $par[1] $par[2] $par[3] $Lconf $par[4] > flush.tmp`;
print `qsub -N $prefix.HCDJ -V -cwd -pe ompi 1 -l virtual_free=32G -o $OUT_calc_HC_JEFF -b y perl $perlPath\'calc_JUNCT_efficiency.pl' $par[0] $par[1] $par[2] $par[3] $conf $par[4] >> flush.tmp`;
print `qsub -N $prefix.HCDJDummy -cwd -V -hold_jid $prefix.HCDJ,$prefix.LCDJ -sync y -b y echo >> flush.tmp`;
print `rm -f flush.tmp`;
print `rm -f $prefix*`;
print "...Done\n";
}




#Generate Intronic Segments
$OUT_INTR_SEGM='Junctions_IntronicSegments.bed';
$OUT_INTR_SEGM_SORTED='Junctions_IntronicSegments_sorted.bed';
unless ($skip{4}){
print "...Generating Intronic Segments bed file\n";
print `perl $perlPath\'tophat_junctions2IntronSegments.pl' Merged_Junctions.bed> $OUT_INTR_SEGM`;
print "...Done\n";
print "...Sorting Intronic Segments File\n";
print `sort -k1,1 -k2,2n $OUT_INTR_SEGM > $OUT_INTR_SEGM_SORTED`;
print "...Done\n";
}




#Calculate Coverage of Intronic Segments by NonJunction Reads:
$NoJunctReads_bed1=$NJoutbed1;		#Generated during Preprocessing
$NoJunctReads_bed2=$NJoutbed2;		#Generated during Preprocessing
$OUT_Coverage_IntrSegm1= $COND1 . '_IntrSegm_coverage.bed';
$OUT_Coverage_IntrSegm2= $COND2 . '_IntrSegm_coverage.bed';
$str_spec="-s";
$str_spec="" if $RNAseq eq "U";
unless ($skip{5}){
print `rm -f *_IntrSegm_coverage.bed`;
print `qsub -N $prefix.IntrCov1 -V -cwd -pe ompi 2 -l virtual_free=32G -o $OUT_Coverage_IntrSegm1 -b y bedtools intersect -c $str_spec -sorted -a $OUT_INTR_SEGM_SORTED -b $NoJunctReads_bed1 > flush.tmp`;
print `qsub -N $prefix.IntrCov2 -V -cwd -pe ompi 2 -l virtual_free=32G -o $OUT_Coverage_IntrSegm2 -b y bedtools intersect -c $str_spec -sorted -a $OUT_INTR_SEGM_SORTED -b $NoJunctReads_bed2 >> flush.tmp`;
print "...Calculating Intronic Segments Read Coverage\n";
print `qsub -N $prefix.IntrCovDummy -cwd -hold_jid $prefix.IntrCov1,$prefix.IntrCov2 -sync y -b y echo >> flush.tmp`;
print `rm -f $prefix.*`;
print `rm -f flush.tmp`;
print "...Done\n";
}




#Calculate Differential Intron Retention  
@par=($OUT_calc_HC_JEFF,$OUT_Coverage_IntrSegm1,$OUT_Coverage_IntrSegm2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc,$Proc_Junctions1,$Proc_Junctions2,"noIRM");
$OUT_INTR_RET="Diff_RetIntr.txt";
unless ($skip{6}){
print `rm -f $OUT_INTR_RET`;
print "...Calculating Differential Intron Retention\n";
print `qsub -N $prefix.DIR -V -cwd -pe ompi 1 -l virtual_free=32G -o $OUT_INTR_RET -b y perl $perlPath\'calc_INTRON_retention.pl' $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] $par[7] $par[8] > flush.tmp`;
print `qsub -N $prefix.DIRDummy -cwd -V -hold_jid $prefix.DIR -sync y -b y echo >> flush.tmp`;
print `rm -f $prefix*`;
print "...Done\n";
}



#Annotate Differential Junctions
@par=($OUT_calc_HC_JEFF,$OUT_calc_LC_JEFF,$olapSel_Junc2Tx,$OUT_INTR_RET,$ENSid2Name,$ENS_Tx_Junc);
$OUT_ANNOT="Annotated_Diff_Junctions.txt";
unless ($skip{7}){
print `rm -f $OUT_ANNOT`;
print "...Annotating Differential Junctions";
print `qsub -N $prefix.ADJ -V -cwd -pe ompi 1 -l virtual_free=32G -o $OUT_ANNOT -b y perl $perlPath\'annotate_Diff_Used_Junctions.pl' $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $conf > flush.tmp`;
print `qsub -N $prefix.ADJDummy -cwd -V -hold_jid $prefix.ADJ -sync y -b y echo >> flush.tmp`;

print `rm -f $prefix*`;
print "...Done\n";
}


if ($IRM eq 'Y'){

#Calculate Differential Intron Retention (IRM mode)
@par=($OUT_calc_LC_JEFF,$OUT_Coverage_IntrSegm1,$OUT_Coverage_IntrSegm2,$olapSel_NJunc12,$olapSel_Junc2Tx,$ENS_Tx_Junc,$Proc_Junctions1,$Proc_Junctions2,"IRM");
$OUT_IRM="Diff_RetIntr_IRM.txt";
unless ($skip{8}){
print `rm -f $OUT_IRM`;
print "...Calculating Differential Intron Retention (IRM mode)\n";
print `qsub -N $prefix.IRM -V -cwd -pe ompi 1 -l virtual_free=32G -o $OUT_IRM -b y perl $perlPath\'calc_INTRON_retention.pl' $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] $par[7] $par[8] > flush.tmp`;
print `qsub -N $prefix.IRMDummy -cwd -V -hold_jid $prefix.IRM -sync y -b y echo >> flush.tmp`;
print `rm -f $prefix*`;
print "...Done\n";
}


#Annotate Differential Introns
@par=($OUT_calc_LC_JEFF,$OUT_calc_LC_JEFF,$olapSel_Junc2Tx,$OUT_IRM,$ENSid2Name,$ENS_Tx_Junc,$SuppJun);
$OUT_IANNOT="Annotated_Diff_Introns.txt";
unless ($skip{9}){
print `rm -f $OUT_IANNOT`;
print "...Annotating Differential Introns";
print `qsub -N $prefix.ADI -V -cwd -pe ompi 1 -l virtual_free=32G -o $OUT_IANNOT -b y perl $perlPath\'annotate_Diff_Used_Introns.pl' $par[0] $par[1] $par[2] $par[3] $par[4] $par[5] $par[6] VHC > flush.tmp`;
print `qsub -N $prefix.ADIDummy -cwd -V -hold_jid $prefix.ADI -sync y -b y echo >> flush.tmp`;

print `rm -f $prefix*`;
print "...Done\n";
}

}


#Build Junction Track










