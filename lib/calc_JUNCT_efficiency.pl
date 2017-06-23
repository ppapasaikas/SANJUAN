use lib $ARGV[6] || "/usr/lib/perl5";
use Text::NSP::Measures::2D::Fisher::twotailed;
use Statistics::Descriptive;


$Proc_Junctions1=$ARGV[0];
$Proc_Junctions2=$ARGV[1];
$olapSel_NJunc12=$ARGV[2];
$olapSel_Junc2Tx=$ARGV[3];
$conf=$ARGV[4];
$ENS_Tx_Junc=$ARGV[5];

$fn_out=$ARGV[7];

$min1read_filter=$ARGV[8];		#Added Claudia 28-10-16
open(OUT,">".$fn_out) || die $!;

open (IN,$Proc_Junctions1)||die;
while (<IN>){
@mat=split /\t/,$_;
$Tot_Jreads_C+=$mat[4];
}
close IN;

open (IN,$Proc_Junctions2)||die;
while (<IN>){
@mat=split /\t/,$_;
$Tot_Jreads_R+=$mat[4];
}
close IN;


$MeanJuncDepth=int(0.5 * ($Tot_Jreads_C+$Tot_Jreads_R));
$NormalRatio1=$MeanJuncDepth/$Tot_Jreads_C;
$NormalRatio2=$MeanJuncDepth/$Tot_Jreads_R;

my %Ccount_unnormalized;

open (IN,$Proc_Junctions1)||die;
while (<IN>){
$line=$_;
@mat=split /\t/,$line;
$JID=$mat[3];
$Ccount{$JID}=int($mat[4]*$NormalRatio1+0.5);	#Normalize Reads to Average Junction Depth
$Ccount_unnormalized{$JID}=$mat[4];
}
close IN;


my %Rcount_unnormalized;

open (IN,$Proc_Junctions2)||die;
while (<IN>){
$line=$_;
@mat=split /\t/,$line;
$JID=$mat[3];
$Rcount{$JID}=int($mat[4]*$NormalRatio2+0.5);	#Normalize Reads to Average Junction Depth
$Rcount_unnormalized{$JID}=$mat[4];
}
close IN;




#Get Neighboring Junctions
open (IN,$olapSel_NJunc12)||die;
while (<IN>){
$line=$_;
@mat=split /\t/,$line;
$JID=$mat[3];
next if $mat[9] eq $JID;
$Tcount{$JID}=$mat[4];
push @{$NEIGH{$JID}},$mat[9];
}
close IN;



#Make Table of All Junct -> Subsuming Txs:
open (IN,$olapSel_Junc2Tx)||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$JID=$mat[3];
next if $mat[9] eq $JID;
push @{$AJn2STxs{$JID}}, $mat[9];	#Junctions 2 Subsuming Txs
$STx{$JID}{$mat[9]}=1;			#Defined entry if the transcript subsumes the JID;
@JFEATS=split("_",$JID);
$JLENG{$JID}=$JFEATS[2]-$JFEATS[1];
$SubJunc{$JID}=1;			#Defined if Junction is subsumed by at least one Transcript;
}
close IN;


open (IN,$ENS_Tx_Junc)||die "\n$ENS_Tx_Junc\n$!\n";
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
push @{$Tx2KnJn{$mat[0]}},$mat[1];
}
close IN;



#Find All immediately adjacent Junctions
foreach $JID (keys %AJn2STxs){
($chr_JID,$start_JID,$end_JID,$strand_JID)=split /\_/,$JID;
	foreach $STx (@{$AJn2STxs{$JID}}){
	$max_start=0;
	$min_end=1e20;
	$UProx="NA";
	$DProx="NA";
		foreach $jn (@{$Tx2KnJn{$STx}}){
		($chr,$start,$end,$strand)=split /\_/,$jn;
		next if ($start==$start_JID && $end==$end_JID);
		if ($start<=$start_JID && $start>$max_start){
			$UProx=$jn;
			$max_start=$start;		
			}		
		if ($end>=$end_JID && $end<$min_end){
			$DProx=$jn;
			$min_end=$end;		
			}
		}
	next if ($DProx eq "NA" && $UProx eq "NA");
    unless ($DProx eq "NA") {$DnProx{$JID}{$DProx}=1;$Prox{$JID}{$DProx}=1}
    unless ($UProx eq "NA") {$UpProx{$JID}{$UProx}=1;$Prox{$JID}{$UProx}=1} 
	#print OUT "$JID\t$UProx\t$DProx\n";	
	}
}



######## Add counts of proper proximal (immediately adjacent) junctions
foreach $JID (keys %Prox){
	foreach $proxj (keys %{$Prox{$JID}} ) {
	next if ($PROPER{$JID}{$proxj});
	next if $JLENG{$proxj}>100000;					#Do not consider spurious/high variance junctions
	next unless ($SubJunc{$proxj});					#Do not consider spurious/high variance junctions
	next if ($Ccount{$proxj}<1 || $Rcount{$proxj}<1);		#Do not consider spurious/high variance junctions
	next if abs(log(($Ccount{$proxj}+1)/($Rcount{$proxj}+1)))>6;	#Do not consider spurious/high variance junctions
	$proper="YES";
		foreach $Tx (@{$AJn2STxs{$proxj}}){	
		$proper="NO" unless $STx{$JID}{$Tx}; #Do not use junction unless it comes from Txs that all of them also subsume JID
		}
		foreach $Tx (@{$AJn2STxs{$JID}}){	
		$proper="NO" unless $STx{$proxj}{$Tx}; #Do not use junction unless it is is subsumed by all JID Tx
		}
	next unless $proper eq "YES";
	$NEIGH_T{$JID}+=$Ccount{$proxj}+$Rcount{$proxj};
	$NEIGH_C{$JID}+=$Ccount{$proxj};
	push @{$C_NEIGH{$JID}},$Ccount{$proxj};
	$NEIGH_R{$JID}+=$Rcount{$proxj};
	push @{$R_NEIGH{$JID}},$Rcount{$proxj};
	$PROPER{$JID}{$proxj}=1;
	}
}



######## Add counts of proper neighboring junctions
foreach $JID (keys %NEIGH){
($chr_JID,$start_JID,$end_JID,$strand_JID)=split /\_/,$JID;
	foreach $neighj ( @{$NEIGH{$JID}} ) {
	($chr,$start,$end,$strand)=split /\_/,$neighj;
	#next if ($start==$start_JID && $end==$end_JID);
	next if $JLENG{$neighj}>100000;					#Do not consider spurious/high variance junctions
	next unless ($SubJunc{$neighj});				#Do not consider spurious/high variance junctions			
	next if abs(log(($Ccount{$neighj}+1)/($Rcount{$neighj}+1)))>6;	#Do not consider spurious/high variance junctions
	$PSI_SNEIGH_R{$JID}+=$Rcount{$neighj} if $start==$start_JID && $end!=$end_JID;
	$PSI_ENEIGH_R{$JID}+=$Rcount{$neighj} if $end==$end_JID && $start!=$start_JID;
	$PSI_SNEIGH_C{$JID}+=$Ccount{$neighj} if $start==$start_JID && $end!=$end_JID;
	$PSI_ENEIGH_C{$JID}+=$Ccount{$neighj} if $end==$end_JID && $start!=$start_JID;
	next if ($Ccount{$neighj}<1 || $Rcount{$neighj}<1);		#Do not consider spurious/high variance junctions
	next if ($PROPER{$JID}{$neighj});
	$proper="YES";
		foreach $Tx (@{$AJn2STxs{$neighj}}){	
		$proper="NO" unless $STx{$JID}{$Tx}; #Do not use junction unless it comes from Txs that all of them also subsume JID
		}
		foreach $Tx (@{$AJn2STxs{$JID}}){	
		$proper="NO" unless $STx{$neighj}{$Tx}; #Do not use junction unless it is is subsumed by all JID Tx
		}
	next unless $proper eq "YES";
	$NEIGH_T{$JID}+=$Ccount{$neighj}+$Rcount{$neighj};
	$NEIGH_C{$JID}+=$Ccount{$neighj};
	push @{$C_NEIGH{$JID}},$Ccount{$neighj};
	$NEIGH_R{$JID}+=$Rcount{$neighj};
	push @{$R_NEIGH{$JID}},$Rcount{$neighj};
	$PROPER{$JID}{$neighj}=1;
	}
}


##### Thresholds: minJNreads, minNGHreads, minJN/NGHreads, maxJNlen, minJNlen, minLFC, Pval, minJNreads per condition)
@VHC= (9,0.100,0.005,100000,50,0.150,0.0001,1);		#Very High Confidence Thresholds
@HC = (7,0.050,0.004,100000,50,0.100,0.0010,1);		#High Confidence Thresholds	(Default)
@MC = (5,0.010,0.002,100000,50,0.050,0.0100,1);		#Medium Confidence Thresholds
@LC = (3,0.001,0.001,200000,50,0.005,0.3000,1);		#Low Confidence Thresholds
@NC = (0,-1e-100,-1e-100,1e100,2,-1e-100,2,0);		#NO Thresholds

@TH=@HC;			# For extremely high-depth data use the VHC settings
@TH=@LC if $conf eq 'LC';
@TH=@MC if $conf eq 'MC';
@TH=@VHC if $conf eq 'VHC';
@TH=@LC if $conf eq 'IRM';
@TH=@NC if $conf eq 'NC';


$,="\t";
$statC = Statistics::Descriptive::Full->new();
$statR = Statistics::Descriptive::Full->new();
foreach $JID (keys %Tcount){
next unless $Tcount{$JID}>$TH[0];			#Sufficient Junction Reads
next unless $NEIGH_T{$JID}>$TH[1] * $Tcount{$JID};	#Sufficient Neighbor Reads

if ($min1read_filter eq 'Y') {		#Added Claudia 28-10-16
	next if ($Ccount{$JID}<$TH[7] || $Rcount{$JID}<$TH[7]);	#Do not consider spurious junctions (minimum 1 junction read in BOTH conditions). 
} elsif ($min1read_filter eq 'N') {
	next if ($Ccount{$JID}<$TH[7] && $Rcount{$JID}<$TH[7]);	# Minimum 1 junction read in AT LEAST ONE conditions (more spurious junctions!!). 
}
next unless ($SubJunc{$JID});				#Do not consider spurious junctions
next unless $Tcount{$JID}>$TH[2] * $NEIGH_T{$JID}; 	#Sufficient Junction/Neighbor Reads (> 0.5%/0.4%/0.1%)
next if $JLENG{$JID}>$TH[3];				#Remove too Long Junctions
next if $JLENG{$JID}<$TH[4];				#Remove too Short Junctions
$Ccount=$Ccount{$JID}+1;
$Rcount=$Rcount{$JID}+1;

$PSI_SNEIGH_C=$PSI_SNEIGH_C{$JID}+1;
$PSI_ENEIGH_C=$PSI_ENEIGH_C{$JID}+1;
$PSI_SNEIGH_R=$PSI_SNEIGH_R{$JID}+1;
$PSI_ENEIGH_R=$PSI_ENEIGH_R{$JID}+1;

$NEIGH_C=$NEIGH_C{$JID}+1;				#Replaced Total Count of Neigh. Junctions with a trimmed mean estimate (03-10-14) for JUNCT_EFF
$NEIGH_R=$NEIGH_R{$JID}+1;
$JE_control=$Ccount/($NEIGH_C+$Ccount);
$JE_transgn=$Rcount/($NEIGH_R+$Rcount);
next unless abs(log(($JE_control/$JE_transgn)))>$TH[5];	 #1.2x/1.1x/1.05x/1.005x fold Threshold -

$statC->clear();
$statR->clear();
$statC->add_data(@{$C_NEIGH{$JID}});
$statR->add_data(@{$R_NEIGH{$JID}});
$tmC = $statC->trimmed_mean(0.4,0.1);
$tmR = $statR->trimmed_mean(0.4,0.1);
#$tmC = $statC->median();
#$tmR = $statR->median();
$NEIGH_TC=int($tmC+1);
$NEIGH_TR=int($tmR+1);

$PSI_S_C=$Ccount/($Ccount+$PSI_SNEIGH_C);
#print OUT "\n$Ccount\t$PSI_SNEIGH_C\t$PSI_ENEIGH_C\n" if $JID=~/13563599/; 
#print OUT "\n$Rcount\t$PSI_SNEIGH_R\t$PSI_ENEIGH_R\n" if $JID=~/13563599/; 
$PSI_E_C=$Ccount/($Ccount+$PSI_ENEIGH_C);
$PSI_S_R=$Rcount/($Rcount+$PSI_SNEIGH_R);
$PSI_E_R=$Rcount/($Rcount+$PSI_ENEIGH_R);
($PSI_S_C,$PSI_S_R)=($PSI_E_C,$PSI_E_R) if $PSI_SNEIGH_C+$PSI_SNEIGH_R<3;
($PSI_E_C,$PSI_E_R)=($PSI_S_C,$PSI_S_R) if $PSI_ENEIGH_C+$PSI_ENEIGH_R<3;
$PSI_C=0.5*($PSI_S_C + $PSI_E_C );
$PSI_R=0.5*($PSI_S_R + $PSI_E_R );


if (($PSI_C-$PSI_R)*($JE_control-$JE_transgn)<0){   ## Original
$PSI_C=$Ccount/$NEIGH_TC;
$PSI_R=$Rcount/$NEIGH_TR;
	if(  ($PSI_C==$PSI_R && $PSI_C==1)|| max($PSI_C,$PSI_R)>1 || ($PSI_C-$PSI_R)*($JE_control-$JE_transgn)<0  ) {
		if ($JE_control>=$JE_transgn){
		$PSI_R=$JE_transgn /$JE_control;
		$PSI_C=1;
		}
		else {
		$PSI_C=$JE_control/$JE_transgn;
		$PSI_R=1;
		}
		if (max($PSI_C,$PSI_R)<1){
		$PSI_C=$PSI_C*0.5*($PSI_C+$PSI_R);
		$PSI_R=$PSI_R*0.5*($PSI_C+$PSI_R);
		}
	}
}

$cj++;

$pval=calculateStatistic(n11=>$Ccount, n1p=>$Ccount+$Rcount, np1=>$NEIGH_C+$Ccount, npp=>$NEIGH_C+$NEIGH_R+$Ccount+$Rcount);
#print OUT "$JID\t$Ccount{$JID}\t$Rcount{$JID}\t$NEIGH_C{$JID}\t$NEIGH_R{$JID}\tTMC: $tmC\tTMR: $tmR\t$pval\t$JLENG{$JID}\t$TH[6]\n"; 

#print OUT "$JID\t$JE_control\t$JE_transgn\t$Ccount{$JID}\t$Rcount{$JID}\t$NEIGH_C{$JID}\t$NEIGH_R{$JID}\t$pval\n" if $pval<$TH[6];
print OUT "$JID\t$JE_control\t$JE_transgn\t$Ccount{$JID}\t$Rcount{$JID}\t$PSI_C\t$PSI_R\t$pval\t$Ccount_unnormalized{$JID}\t$Rcount_unnormalized{$JID}\n" if $pval<$TH[6];
}
print OUT ">$cj\n";

close(OUT);




# wait 5 min so that output on cluster is written completely for sure 
sleep 60;





sub min{
@vals=@_;
$min=1e20;
foreach (@vals){
$min=$_ if $_<$min;
}
return $min;
}

sub max{
@vals=@_;
$max=-1e20;
foreach (@vals){
$max=$_ if $_>$max; 
}
return $max;
}










