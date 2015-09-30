$file=$ARGV[0];
$helpful_message="\nUsage: perl SANJUAN_analyze_CEx.pl SANJUAN_RESULTS_FILE Name_of_Condition1 Name_of_Condition2\n";
if ($#ARGV<0){
die "$helpful_message";
}

####### Open File and get condition names ########


if ($#ARGV<2){
die "$helpful_message";
}
else {
$cond1=$ARGV[1];
$cond2=$ARGV[2];
}

$JuncFile1="Junctions_" . $cond1 . '.bed';
$JuncFile2="Junctions_" . $cond2 . '.bed';

open (JF1,"$JuncFile1") || die $! . "\n$Juncfile1\n";
while (<JF1>) {
$l=$_;
chomp $l;
@m=split /\t/,$l;
$JSTARTS{$m[3]}=$m[1];
$JENDS{$m[3]}=$m[2];
$JSTARTING{$m[1]}{$m[3]}++;
$JENDING{$m[2]}{$m[3]}++;
$JCOUNT1{$m[3]}=$m[4];
}

open (JF2,"$JuncFile2") || die $! . "\n$Juncfile2\n";
while (<JF2>) {
$l=$_;
chomp $l;
@m=split /\t/,$l;
$JSTARTS{$m[3]}=$m[1];
$JENDS{$m[3]}=$m[2];
$JSTARTING{$m[1]}{$m[3]}++;
$JENDING{$m[2]}{$m[3]}++;
$JCOUNT2{$m[3]}=$m[4];
$Delta{$m[3]}=$JCOUNT2{$m[3]}-$JCOUNT1{$m[3]};
}

### %JSTARTS %JENDS %JSTARTING %JENDING %Delta






###################### Define differential JUNCTIONS, GENE and corresponding PSIs ######################
open (FILE,"$file") || die $! . "\n$file\n$helpful_message\n";
while (<FILE>) {
$l=$_;
chomp $l;
@m=split /\t/,$l;

$gene=$m[1];			
$hcj=$m[2];			
$compj=$m[3];
$psi1=$m[22];
$psi2=$m[23];
$cpsi1=$m[28];
$cpsi2=$m[29];
next unless $hcj=~ /chr/;    #To remove headers

###################### Define COMPETITION TYPE (if needed) ######################
$comptype=$m[5];
#if ($comptype=~ /DUAL_COMP/) {$type="DUAL_COMP"};
##else if ($comptype=~ /INTRON/) {$type="RET_INTRON"};
#else if ($comptype=~ /CE/) {$type="CE"};
#else if ($comptype=~ /5'SS/) {$type="5'SS"};
#else if ($comptype=~ /3'SS/) {$type="3'SS"};

###################### Define LENGTHS of the JUNCTIONS ######################
@coordhcj=split /_/, $hcj;
@coordcompj=split /_/, $compj;     #split both junctions ids by _
$lenhcj=$coordhcj[2]-$coordhcj[1];        
$lencompj=$coordcompj[2]-$coordcompj[1];    #length is 3rd-2nd split element


###################### Define relationships ######################

$LEN{$hcj}=$lenhcj;		#junctions to their length
$LEN{$compj}=$lencompj;

$JUNC2GENE{$hcj}=$gene;		#junctions to gene name
$JUNC2GENE{$compj}=$gene;

$JUNC2PSI{$hcj}{$cond1}=$psi1;		#junctions to PSI in a particular stage 
$JUNC2PSI{$hcj}{$cond2}=$psi2;		#if you use these PSIs, you will only have the PSIs from the significant junctions in each comparison 
$JUNC2PSI{$compj}{$cond1}=$cpsi1;
$JUNC2PSI{$compj}{$cond2}=$cpsi2;

$JUNCPAIR2TYPE{$hcj}{$compj}=$comptype;		#junction pairs (HCJ-competing) to competition type
$JUNCPAIR2TYPE{$compj}{$hcj}=$comptype;
}

close FILE;



print "SKIP_JUNC\tCHROMOSOME\tE_START\tE_END\tSTRAND\tGENE_SYMBOL\tPSI_COND1\tPSI_COND2\tDELTA_PSI_(COND2-COND1)\n";

foreach $jid (keys %JUNC2GENE) {				#foreach junctionid (contained in JUNC2GENE) {
next if $saw{$jid}>0;							#only unique junction values (if junction was already seen, go to the next step)
$saw{$jid}++;
  foreach $pairj (keys %{$JUNCPAIR2TYPE{$jid}}  ) {		#foreach of their competing junctions {
  $c=0;
  next if $ssaw{$jid}{$pairj}>0; 
  $ssaw{$jid}{$pairj}++;
  $ssaw{$pairj}{$jid}++;				
  next unless $JUNCPAIR2TYPE{$jid}{$pairj}=~/CE/;		#check if junction pair is of ***type CE***
  $minDeltaDist=1e20;
  $exon_start="NA";
  $exon_end="NA";
	if ($LEN{$jid}<$LEN{$pairj}) {				#check which pair element is the shortest and 
	@coordjid=split /_/, $jid;	
	@coordpairj=split /_/, $pairj;
	$chr=$coordjid[0];
	$strand=$coordjid[3];
		if ($coordpairj[1]<$coordjid[1]){
		$exon_end=$coordjid[1]-1;
		$target_start=$coordpairj[1];
		$target_max=$coordjid[1]-1;
		$target_delta=$CountCond2{$jid}-$CountCond1{$jid};		
			foreach $tjid (keys %{$JSTARTING{$target_start}}){
			next unless $JENDS{$tjid}<$target_max;			
			$DeltaDist=$Delta{$tjid}-$Delta{$jid};
			if ($DeltaDist<$minDeltaDist) {$minDeltaDist=$DeltaDist;$TargetJunction=$tjid;@coordtjid=split /_/,$tjid;$exon_start=$coordtjid[2]+1;}
			}		
		}

		else {
		$exon_start=$coordjid[2]+1;
		$target_end=$coordpairj[2];
		$target_min=$coordjid[2]+1;
		$target_delta=$CountCond2{$jid}-$CountCond1{$jid};		
			foreach $tjid (keys %{$JENDING{$target_end}}){
			next unless $JSTARTS{$tjid}>$target_min;			
			$DeltaDist=$Delta{$tjid}-$Delta{$jid};
			if ($DeltaDist<$minDeltaDist) {$minDeltaDist=$DeltaDist;$TargetJunction=$tjid;@coordtjid=split /_/,$tjid;$exon_end=$coordtjid[1]-1;}
			}		
		}

	$DPSI=$JUNC2PSI{$jid}{$cond2}-$JUNC2PSI{$jid}{$cond1};
	print "$pairj\t$chr\t$exon_start\t$exon_end\t$strand\t$JUNC2GENE{$jid}\t\t$JUNC2PSI{$jid}{$cond1}\t$JUNC2PSI{$jid}{$cond2}\t$DPSI";	#print PSIs of shortest across all the stages
	}





	else  {
	@coordjid=split /_/, $jid;	
	@coordpairj=split /_/, $pairj;
	$chr=$coordjid[0];
	$strand=$coordjid[3];
		if ($coordjid[1]<$coordpairj[1]){
		$exon_end=$coordpairj[1]-1;
		$target_start=$coordjid[1];
		$target_max=$coordpairj[1]-1;
		$target_delta=$CountCond2{$pairj}-$CountCond1{$pairj};		
			foreach $tjid (keys %{$JSTARTING{$target_start}}){
			next unless $JENDS{$tjid}<$target_max;			
			$DeltaDist=$Delta{$tjid}-$Delta{$pairj};
			if ($DeltaDist<$minDeltaDist) {$minDeltaDist=$DeltaDist;$TargetJunction=$tjid;@coordtjid=split /_/,$tjid;$exon_start=$coordtjid[2]+1;}
			}		
		}

		else {
		$exon_start=$coordpairj[2]+1;
		$target_end=$coordjid[2];
		$target_min=$coordpairj[2]+1;
		$target_delta=$CountCond2{$pairj}-$CountCond1{$pairj};		
			foreach $tjid (keys %{$JENDING{$target_end}}){
			next unless $JSTARTS{$tjid}>$target_min;			
			$DeltaDist=$Delta{$tjid}-$Delta{$pairj};
			if ($DeltaDist<$minDeltaDist) {$minDeltaDist=$DeltaDist;$TargetJunction=$tjid;@coordtjid=split /_/,$tjid;$exon_end=$coordtjid[1]-1;}
			}		
		}

	$DPSI=$JUNC2PSI{$pairj}{$cond2}-$JUNC2PSI{$pairj}{$cond1};
	print "$jid\t$chr\t$exon_start\t$exon_end\t$strand\t$JUNC2GENE{$pairj}\t$JUNC2PSI{$pairj}{$cond1}\t$JUNC2PSI{$pairj}{$cond2}\t$DPSI";
	}
  print "\n";
  }
}






########################################################################################    

sub median {							#define median function
my ($ar) = shift;
my $count = scalar @$ar;
# Sort a COPY of the array, leaving the original untouched
my @array = sort { $a <=> $b } @$ar;
if ($count % 2) {
return $array[int($count/2)];
} else {
return ($array[$count/2] + $array[$count/2 - 1]) / 2;
}
}








