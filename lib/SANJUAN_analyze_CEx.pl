$helpful_message="\nUsage: perl SANJUAN_analyze_CEx.pl SANJUAN_RESULTS_FILE Diff_Junctions_NC.txt\n";


####### Open File and get condition names ########

if ($#ARGV<1){
die "$helpful_message";
}
else {
$RESfile=$ARGV[0];
$NCfile=$ARGV[1];
}




###################### Define differential JUNCTIONS, GENE and corresponding PSIs ######################
open (FILE,"$RESfile") || die $! . "\n$file\n$helpful_message\n";
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

if ($l=~/^INCL/){
$title=$l;
@mat=split /\t/,$title;
@matt=@mat;
$cond1=$1 if $mat[18]=~/Eff_(.+)/;
$cond2=$1 if $mat[19]=~/Eff_(.+)/;
}



next unless $hcj=~ /chr/;    #To remove headers

###################### Define COMPETITION TYPE (if needed) ######################
$comptype=$m[5];
#if ($comptype=~ /DUAL_COMP/) {$type="DUAL_COMP"};
##else if ($comptype=~ /INTRON/) {$type="RET_INTRON"};
#else if ($comptype=~ /CE/) {$type="CE"};
#else if ($comptype=~ /5'SS/) {$type="5'SS"};
#else if ($comptype=~ /3'SS/) {$type="3'SS"};
next unless $comptype=~/CE/;		#check if junction pair is of ***type CE***

###################### Define LENGTHS of the JUNCTIONS ######################
@coordhcj=split /_/, $hcj;
@coordcompj=split /_/, $compj;     #split both junctions ids by _
$lenhcj=$coordhcj[2]-$coordhcj[1];        
$lencompj=$coordcompj[2]-$coordcompj[1];    #length is 3rd-2nd split element


###################### Define relationships ######################

$LEN{$hcj}=$lenhcj;		#junctions to their length
$LEN{$compj}=$lencompj;

if ($LEN{$hcj}>$LEN{$compj}){
$LONGERJ{$hcj}=1;
}
else {
$LONGERJ{$compj}=1;
}


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






###################### Define differential JUNCTIONS, GENE and corresponding PSIs ######################
open (FILE,"$NCfile") || die $! . "\n$file\n$helpful_message\n";
while (<FILE>) {
$l=$_;
next unless $l=~/^chr/;	#Remove headers
chomp $l;
@m=split /\t/,$l;

$jid=$m[0];			
$jeff1=$m[1];			
$jeff2=$m[2];
$psi1=$m[5];
$psi2=$m[6];
@coordjid=split /_/, $jid;
$jstart=$coord[0].$coordjid[1].$coord[3];
$jend=$coord[0].$coordjid[2].$coord[3];
$JSTARTS{$jstart}{$jid}=1;
$JENDS{$jend}{$jid}=1;


$NCJUNC2PSI1{$jid}=sprintf("%.3f", $psi1);		#junctions to PSI in a particular stage 
$NCJUNC2PSI2{$jid}=sprintf("%.3f",$psi2);
$NCJUNC2JEFF1{$jid}=sprintf("%.3f",$jeff1);
$NCJUNC2JEFF2{$jid}=sprintf("%.3f",$jeff2);


}

close FILE;





print "SKIP_JUNC\tCHROMOSOME\tE_START\tE_END\tSTRAND\tGENE_SYMBOL\tPSI_$cond1\tPSI_$cond2\tDELTA_PSI_($cond2-$cond1)\tJEFF_$cond1\tJEFF_$cond2\tDELTA_JEFF_($cond2-$cond1)\n";


foreach $jid (keys %LONGERJ) {				#foreach junctionid (contained in JUNC2GENE) {
next if $saw{$jid}>0;						#only unique junction values (if junction was already seen, go to the next step)
$saw{$jid}++;
@coordjid=split /_/, $jid;	
$jstart=$coord[0].$coordjid[1].$coord[3];
$jend=$coord[0].$coordjid[2].$coord[3];
$minDeltaDist=1e20;
$DELTAPSI_LJ=$NCJUNC2PSI2{$jid}-$NCJUNC2PSI1{$jid};

  foreach $pairj (keys   %{$JSTARTS{$jstart}} ) {		#foreach of their competing junctions {
  $c=0;
  next if $pairj eq $jid;
  $DELTAPSI_PJ=$NCJUNC2PSI2{$pairj}-$NCJUNC2PSI1{$pairj};			
  $DeltaDist=abs($DELTAPSI_LJ+$DELTAPSI_PJ);
  $DeltaDist=2*(0.01+$DeltaDist) if $DELTAPSI_PJ*$DELTAPSI_LJ>=0 || abs($DELTAPSI_PJ)<0.5*abs($DELTAPSI_LJ);
	#print "$jid\t$pairj\tDPSIL:$DELTAPSI_LJ\t$DELTAPSI_PJ\n";
  	if ($DeltaDist<$minDeltaDist){
  	$Spairj=$pairj;
  	$minDeltaDist=$DeltaDist;
  	}
  }
next if $minDeltaDist>10;


$minDeltaDist=1e20;

  foreach $pairj (keys   %{$JENDS{$jend}} ) {		#foreach of their competing junctions {
  $c=0;
  next if $pairj eq $jid;
  $DELTAPSI_PJ=$NCJUNC2PSI2{$pairj}-$NCJUNC2PSI1{$pairj};			
  $DeltaDist=abs($DELTAPSI_LJ+$DELTAPSI_PJ);
  $DeltaDist=2*(0.01+$DeltaDist) if $DELTAPSI_PJ*$DELTAPSI_LJ>=0 || abs($DELTAPSI_PJ)<0.5*abs($DELTAPSI_LJ);		
  #$DeltaDist=abs($NCJUNC2PSI2{$pairj}-$NCJUNC2PSI1{$pairj}+$NCJUNC2PSI2{$jid}-$NCJUNC2PSI1{$jid});
	#print "$jid\t$pairj\tDPSIL:$DELTAPSI_LJ\t$DELTAPSI_PJ\n";
  	if ($DeltaDist<$minDeltaDist){
  	$Epairj=$pairj;
 	$minDeltaDist=$DeltaDist;
  	}
  }
next if $minDeltaDist>10;
	
@coordjid=split /_/, $jid;	
@coordSpairj=split /_/, $Spairj;
@coordEpairj=split /_/, $Epairj;
$chr=$coordjid[0];
$strand=$coordjid[3];



$exon_start=$coordSpairj[2]+1;
$exon_end=$coordEpairj[1]-1;
next if $exon_start>=$exon_end;
$PSI1=0.5*$NCJUNC2PSI1{$Spairj}+0.5*$NCJUNC2PSI1{$Epairj};	#Correct this (when do we need to substract or not, plus printing of JEFF values)
$PSI2=0.5*$NCJUNC2PSI2{$Spairj}+0.5*$NCJUNC2PSI2{$Epairj};
$JEFF1=0.5*$NCJUNC2JEFF1{$Spairj}+0.5*$NCJUNC2JEFF1{$Epairj};
$JEFF2=0.5*$NCJUNC2JEFF2{$Spairj}+0.5*$NCJUNC2JEFF2{$Epairj};
$DPSI=sprintf("%.3f",$PSI2-$PSI1);
$DJEFF=sprintf("%.3f",$JEFF2-$JEFF1);
print "$jid\t$chr\t$exon_start\t$exon_end\t$strand\t$JUNC2GENE{$jid}\t$PSI1\t$PSI2\t$DPSI\t$JEFF1\t$JEFF2\t$DJEFF\n";	#print PSIs of shortest across all the stages


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








