use lib $ARGV[9] || "/usr/lib/perl5";
use Text::NSP::Measures::2D::Fisher::twotailed;


$Diff_Junct_Eff=$ARGV[0];
$Coverage_IntrSegm1=$ARGV[1];
$Coverage_IntrSegm2=$ARGV[2];
$olapSel_NJunc12=$ARGV[3];
$olapSel_Junc2Tx=$ARGV[4];
$ENS_Tx_Junc=$ARGV[5];
$Proc_Junctions1=$ARGV[6];
$Proc_Junctions2=$ARGV[7];
$IRM=$ARGV[8];

open (IN, $Diff_Junct_Eff) ||die;
while (<IN>){
$line=$_;
next if $line=~/^>/;
chomp $line;
@mat=split /\t/,$line;
$LR{$mat[0]}=log(($mat[2]/$mat[1]));
$DiffJunct{$mat[0]}=defined;
}
close IN;


if ($IRM ne 'IRM'){
open (IN, $Coverage_IntrSegm1)||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$Tot_Ireads_C+=$mat[6];
}
close IN;

open (IN, $Coverage_IntrSegm2)||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$Tot_Ireads_R+=$mat[6];
}
close IN;

$MeanIntrSegDepth=int(0.5 * ($Tot_Ireads_C+$Tot_Ireads_R));
$NormalRatio1=$MeanIntrSegDepth/$Tot_Ireads_C;
$NormalRatio2=$MeanIntrSegDepth/$Tot_Ireads_R;
}

else {
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
}



open (IN, $Coverage_IntrSegm1)||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$SID=$mat[3];
$IID=$1 if $SID=~/(.+[\-\+])/;
$S_Ccount{$SID}=int($mat[6]*$NormalRatio1+0.5);		#Normalize Reads to Control Intron Coverage
$I_Ccount{$IID}+=int($mat[6]*$NormalRatio1+0.5);	#Normalize Reads to Control Intron Coverage
$S_Tcount{$SID}=$S_Ccount{$SID};
$I_Tcount{$IID}+=$I_Ccount{$IID};
$Nsegm{$IID}++;
}
close IN;


open (IN, $Coverage_IntrSegm2)||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$SID=$mat[3];
$IID=$1 if $SID=~/(.+[\-\+])/;
$S_Rcount{$SID}=int($mat[6]*$NormalRatio2+0.5);		#Normalize Reads to Control Intron Coverage
$I_Rcount{$IID}+=int($mat[6]*$NormalRatio2+0.5);	#Normalize Reads to Control Intron Coverage
$S_Tcount{$SID}+=$S_Rcount{$SID};
$I_Tcount{$IID}+=$I_Rcount{$IID};
}
close IN;


open (IN,$Proc_Junctions1)||die;
while (<IN>){
$line=$_;
@mat=split /\t/,$line;
$JID=$mat[3];
$Ccount{$JID}=int($mat[4]*$NormalRatio1+0.5);	#Normalize Reads to Average Junction Depth
}
close IN;

open (IN,$Proc_Junctions2)||die;
while (<IN>){
$line=$_;
@mat=split /\t/,$line;
$JID=$mat[3];
$Rcount{$JID}=int($mat[4]*$NormalRatio2+0.5);	#Normalize Reads to Average Junction Depth
}
close IN;





#Get Neighboring Junctions
open (IN,$olapSel_NJunc12)||die;
while (<IN>){
$line=$_;
@mat=split /\t/,$line;
$JID=$mat[3];
next if $mat[9] eq $JID;
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


open (IN,$ENS_Tx_Junc)||die;
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
	#print "$JID\t$UProx\t$DProx\n";	
	}
}



######## Add counts of proper proximal (immediately adjacent) junctions (here -> introns)
foreach $JID (keys %Prox){
	foreach $proxj (keys %{$Prox{$JID}} ) {
	next if ($PROPER{$JID}{$proxj});
	next if $JLENG{$proxj}>100000;					#Do not consider spurious/high variance junctions
	next unless ($SubJunc{$proxj});					#Do not consider spurious/high variance junctions
	next if ($I_Ccount{$proxj}<1 || $I_Rcount{$proxj}<1);		#Do not consider spurious/high variance junctions
	next if abs(log(($I_Ccount{$sproxj}+1)/($I_Rcount{$proxj}+1)))>6;	#Do not consider spurious/high variance junctions
	$proper="YES";
		foreach $Tx (@{$AJn2STxs{$proxj}}){	
		$proper="NO" unless $STx{$JID}{$Tx}; #Do not use junction unless it comes from Txs that all of them also subsume JID
		}
		foreach $Tx (@{$AJn2STxs{$JID}}){	
		$proper="NO" unless $STx{$proxj}{$Tx}; #Do not use junction unless it is is subsumed by all JID Tx
		}
	next unless $proper eq "YES";
	if ($IRM ne 'IRM'){
	$NEIGH_C{$JID}+=$I_Ccount{$proxj};
	$NEIGH_R{$JID}+=$I_Rcount{$proxj};
	}
	else{
	$NEIGH_C{$JID}+=$Ccount{$proxj};
	$NEIGH_R{$JID}+=$Rcount{$proxj};
	}
	$PROPER{$JID}{$proxj}=1;
	}
}

######## Add counts of proper neighboring junctions (here -> introns)
foreach $JID (keys %NEIGH){
	foreach $neighj ( @{$NEIGH{$JID}} ) {
	next if ($PROPER{$JID}{$neighj});
	next if $JLENG{$neighj}>100000;					#Do not consider spurious/high variance junctions
	next unless ($SubJunc{$neighj});				#Do not consider spurious/high variance junctions			
	next if ($I_Ccount{$neighj}<1 || $I_Rcount{$neighj}<1);		#Do not consider spurious/high variance junctions
	next if abs(log(($I_Ccount{$neighj}+1)/($I_Rcount{$neighj}+1)))>6;	#Do not consider spurious/high variance junctions
	$proper="YES";
		foreach $Tx (@{$AJn2STxs{$neighj}}){	
		$proper="NO" unless $STx{$JID}{$Tx}; #Do not use junction unless it comes from Txs that all of them also subsume JID
		}
		foreach $Tx (@{$AJn2STxs{$JID}}){	
		$proper="NO" unless $STx{$neighj}{$Tx}; #Do not use junction unless it is is subsumed by all JID Tx
		}
	next unless $proper eq "YES";
	if ($IRM ne 'IRM'){
	$NEIGH_C{$JID}+=$I_Ccount{$neighj};
	$NEIGH_R{$JID}+=$I_Rcount{$neighj};
	}
	else {	
	$NEIGH_C{$JID}+=$Ccount{$neighj};
	$NEIGH_R{$JID}+=$Rcount{$neighj};
	}	
	$PROPER{$JID}{$neighj}=1;
	}
}





@VHC=(9,0.200,0.005,100000,50,0.180,0.000001);		#Very High Confidence Thresholds
@HC =(7,0.050,0.004,100000,50,0.100,0.001000);		#High Confidence Thresholds
@LC =(3,0.001,0.001,200000,50,0.005,0.300000);		#Low Confidence Thresholds
@TH=@HC;
$pvalTH=0.01;

foreach $SID (keys %S_Tcount){
$IID=$1 if $SID=~/(.+[\-\+])/;
next unless $DiffJunct{$IID} || $IRM eq 'IRM';
next if ($S_Ccount{$SID}+$S_Rcount) < 0.05 * $I_Tcount{$IID}/$Nsegm{$IID} && ($S_Ccount{$SID}+$S_Rcount)< 5; #Discard if segment has low coverage
$pval=1;
$Delta=0;

$S_Ccount=$S_Ccount{$SID}+1;
$S_Rcount=$S_Rcount{$SID}+1;
$NEIGH_C=$NEIGH_C{$IID}+1;
$NEIGH_R=$NEIGH_R{$IID}+1;
$SR_control=$S_Ccount/($NEIGH_C+$S_Ccount);
$SR_transgn=$S_Rcount/($NEIGH_R+$S_Rcount);
$IRLR=log($SR_transgn/$SR_control);
$pval=calculateStatistic(n11=>$S_Ccount, n1p=>$S_Ccount+$S_Rcount, np1=>$NEIGH_C+$S_Ccount, npp=>$NEIGH_C+$NEIGH_R+$S_Ccount+$S_Rcount);

#print "$IID\t$S_Ccount\t$S_Rcount\t$NEIGH_C\t$NEIGH_R\t$IRLR\t$pval\n" if $IID=~/18168624/;

push @{$IRLRS{$IID}},$IRLR;
push @{$PVALS{$IID}},$pval;
}


$,="\t";
foreach $IID (keys %IRLRS){
$medIRLR=median(\@{$IRLRS{$IID}});
$medPval=geomPval(\@{$PVALS{$IID}},\@{$IRLRS{$IID}});	#Changed on 03-2015
print "$IID\t$I_Ccount{$IID}\t$I_Rcount{$IID}\t$NEIGH_C{$IID}\t$NEIGH_R{$IID}\t$LR{$IID}\t$medIRLR\t$medPval\n";# if $medPval<$pvalTH;
}




# wait 5 min so that output on cluster is written completely for sure 
sleep 300;



sub median {
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


sub medianPval {	#Median Pvalue function
my ($ar1) = $_[0];
my ($ar2) = $_[1];
@pvals=@$ar1;
@vals=@$ar2;
$medVal=median (\@vals);
	for $c (0..$#pvals){
	$pvals[$c]=0.5 unless $medVal*$vals[$c]>0;
	}
return median(\@pvals);
}


sub geomPval {		#Geometric mean Pvalue function
my ($ar1) = $_[0];
my ($ar2) = $_[1];
my $LP=0;
@pvals=@$ar1;
@vals=@$ar2;
$medVal=median (\@vals);
	for $c (0..$#pvals)
	{	
	$pvals[$c]=1 unless $medVal*$vals[$c]>0;
	}

	foreach (@pvals){ 
	#next if $_ < 1e-100;	#bug corrected 04-15:
	if ($_ < 1e-100){$add=-230}	#ln(10^-100) ~ -230
	else {$add=log($_)}
	$LP=$LP + $add; 
	}
$gmP=exp($LP/($#pvals+1));
return $gmP;
}









