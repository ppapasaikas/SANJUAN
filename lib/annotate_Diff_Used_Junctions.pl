$Diff_Junct_Eff_HC=$ARGV[0];
$Diff_Junct_Eff_LC=$ARGV[1];
$olapSel_Junc2Tx=$ARGV[2];
$Diff_Intr_Ret=$ARGV[3];
$ENSid2Name=$ARGV[4];
$ENS_Tx_Junc=$ARGV[5];
$conf=$ARGV[6];
$COND1=$ARGV[7];
$COND2=$ARGV[8];
$fn_out=$ARGV[9];
open(OUT,">".$fn_out) || die $!;



($minDPSI,$minPvRET,$minLFC)=(0.20,0.01,0.4) if $conf eq 'VHC'; #Pan: !!Note!! Only the minDPSI is relevant for junctions. minPvRET and minLFC are only applied to intron retention (Non-IRM mode)
($minDPSI,$minPvRET,$minLFC)=(0.15,0.05,0.2) if $conf eq 'HC';
($minDPSI,$minPvRET,$minLFC)=(0.10,0.10,0.1) if $conf eq 'MC';
($minDPSI,$minPvRET,$minLFC)=(0,0.3,0.005) if $conf eq 'LC';  ## added Andre Jun 21, 2017 for getting LC junctions and later for getting LC CEs
($minDPSI,$minPvRET,$minLFC)=(0,2,0) if $conf eq 'NC'; ## added Andre Sep 8, 2016 for getting all exons / introns


open (IN, $ENSid2Name)||die "\n$ENSid2Name\n$!\n";
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$ENSid2Gname{$mat[0]}=$mat[1];
}
close IN;


open (IN,$ENS_Tx_Junc)||die "\n$!\n";
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
next if ($ENS_SAW{$mat[1]});
$ENS_SAW{$mat[1]}=defined;
@mm=split /\_/,$mat[1];
$CHR=$mm[0];
$STRAND=$mm[3];
$STARTJ=$CHR . '_' . $mm[1] . '_' . $STRAND;
$ENDJ=$CHR . '_' . $mm[2] . '_' . $STRAND;
$JUNCTION{$mat[1]}++;
$DONOR{$STARTJ}++ if $STRAND eq '+';
$DONOR{$ENDJ}++ if $STRAND eq '-';
$ACCPT{$ENDJ}++ if $STRAND eq '+';
$ACCPT{$STARTJ}++ if $STRAND eq '-';

$EL_CHR{$mat[1]}=$mm[0];
$EL_START{$mat[1]}=$mm[1];
$EL_END{$mat[1]}=$mm[2];
$EL_STRAND{$mat[1]}=$mm[3];
$seg=1+int(($mm[1]+$mm[2])/200000);
$loc=$mm[0] . $mm[3] . $seg;
push @{$LOCJ{$loc}},$mat[1];
}
close IN;


open (IN, $Diff_Junct_Eff_HC) ||die;
while (<IN>){
$line=$_;
next if $line=~/^>/;
chomp $line;
@mat=split /\t/,$line;
@mm=split /\_/,$mat[0];
next unless abs($mat[6]-$mat[5])>=$minDPSI; ## changed from > to >= Andre Sep 8 for getting all exons / introns	
@{$ATTRIBS{$mat[0]}}=@mat;
$DELTA{$mat[0]}=$mat[6]-$mat[5];
$LR{$mat[0]}=log(($mat[2]/$mat[1]));
$CHR{$mat[0]}=$mm[0];
$START{$mat[0]}=$mm[1];
$END{$mat[0]}=$mm[2];
$STRAND{$mat[0]}=$mm[3];
}
close IN;




open (IN, $Diff_Junct_Eff_LC) ||die;
while (<IN>){
$line=$_;
next if $line=~/^>/;
chomp $line;
@mat=split /\t/,$line;
@mm=split /\_/,$mat[0];
@{$L_ATTRIBS{$mat[0]}}=@mat;
$L_LR{$mat[0]}=log(($mat[2]/$mat[1]));
$L_DELTA{$mat[0]}=$mat[6]-$mat[5];
$EL_CHR{$mat[0]}=$mm[0];
$EL_START{$mat[0]}=$mm[1];
$EL_END{$mat[0]}=$mm[2];
$EL_STRAND{$mat[0]}=$mm[3];
$seg=1+int(($mm[1]+$mm[2])/200000);
$loc=$mm[0] . $mm[3] . $seg;
push @{$LOCJ{$loc}},$mat[0];
}
close IN;


open (IN, $olapSel_Junc2Tx)||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
next unless $ATTRIBS{$mat[3]};

	unless ($ENSid2Gname{$mat[9]}){
	${$Gnames{$mat[3]}}[0]="NA";
	next;
	}

	if (${$Gnames{$mat[3]}}[0] eq "NA") {${$Gnames{$mat[3]}}[0]=$ENSid2Gname{$mat[9]}} else {
	push @{$Gnames{$mat[3]}},$ENSid2Gname{$mat[9]}
	}
}
close IN;






open (IN, $Diff_Intr_Ret) ||die;
while (<IN>){
$line=$_;
next if $line=~/^>/;
chomp $line;
@mat=split /\t/,$line;
next unless ($mat[5]*$mat[6])<0;	#Supporting Junction Evidence
next unless $mat[7] < $minPvRET;	#below Pval threshold
next if abs($mat[6])<$minLFC;		#above logfold threshold
$RET_INTR{$mat[0]}=defined;
@{$RI{$mat[0]}}=@mat;
}
close IN;




#Junctions with same start/end using the UCSC Ensembl Dataset:
foreach $ej1 (keys %JUNCTION){
$seg1=0+int(($EL_START{$ej1}+$EL_END{$ej1})/200000);
$seg2=$seg1+1;
$seg3=$seg1+2;
$loc1=$EL_CHR{$ej1} . $EL_STRAND{$ej1} . $seg1;
$loc2=$EL_CHR{$ej1} . $EL_STRAND{$ej1} . $seg2;
$loc3=$EL_CHR{$ej1} . $EL_STRAND{$ej1} . $seg3;
	foreach $ej2 (@{$LOCJ{$loc1}},@{$LOCJ{$loc2}},@{$LOCJ{$loc3}}){
	next if $ej1 eq $ej2;
	push @{$Same_Start{$ej1}},$ej2 if $EL_START{$ej1}==$EL_START{$ej2};
	push @{$Same_End{$ej1}},$ej2 if $EL_END{$ej1}==$EL_END{$ej2};
	}
undef %saw;
@{$Same_Start{$ej1}} = grep(!$saw{$_}++, @{$Same_Start{$ej1}});
@{$Same_End{$ej1}} = grep(!$saw{$_}++, @{$Same_End{$ej1}});
}

#Junctions with same start/end using the Low Confidence Dataset:
foreach $lcj1 (keys %L_ATTRIBS){
$seg1=0+int(($EL_START{$lcj1}+$EL_END{$lcj1})/200000);
$seg2=$seg1+1;
$seg3=$seg1+2;
$loc1=$EL_CHR{$lcj1} . $EL_STRAND{$lcj1} . $seg1;
$loc2=$EL_CHR{$lcj1} . $EL_STRAND{$lcj1} . $seg2;
$loc3=$EL_CHR{$lcj1} . $EL_STRAND{$lcj1} . $seg3;
	foreach $lcj2 (@{$LOCJ{$loc1}},@{$LOCJ{$loc2}},@{$LOCJ{$loc3}}){
	next if $lcj1 eq $lcj2;
	push @{$Same_Start{$lcj1}},$lcj2 if $EL_START{$lcj1}==$EL_START{$lcj2};
	push @{$Same_End{$lcj1}},$lcj2 if $EL_END{$lcj1}==$EL_END{$lcj2};
	}
undef %saw;
@{$Same_Start{$lcj1}} = grep(!$saw{$_}++, @{$Same_Start{$lcj1}});
@{$Same_End{$lcj1}} = grep(!$saw{$_}++, @{$Same_End{$lcj1}});
}




$PROXIMAL{"NA"}="NA";
$L_LR{"NA"}="NA";
$L_ATTRIBS{"NA"}="NA";


# We search for pairs of competing junctions 
# - efficiency of one goes up and efficiency of the other down to a very similar degree
# - both junctions are close (each junction is a pair if start-,end-coordinate, (minimal) distance 
#   between 2 junctions is the smallest distance between each combination of their start/end coordinates
# - we start with all high confidence junctions and search for a competing junction among all low confidence junctions
# For each identified pair of junction and its competing junction we determine the type of alternative splicing event 

foreach $hcj (keys %ATTRIBS){   # go over all high conf. junctions
$Dist{$hcj}=7.0;
$TYPE{$hcj}="UNIDNT_COMP";
$PROXIMAL{$hcj}="NA";
$PROX_DIS{$hcj}="NA";
$LRdist{$hcj}="NA";
$PVdist{$hcj}="NA";
$RCdist{$hcj}="NA";
$DSdist{$hcj}="NA";
$PROX_ABS_DIS{$hcj}=$minDist;
$seg1=0+int(($START{$hcj}+$END{$hcj})/200000);
$seg2=$seg1+1;
$seg3=$seg1+2;
$loc1=$CHR{$hcj} . '-' . $seg1;
$loc2=$CHR{$hcj} . '-' . $seg2;
$loc3=$CHR{$hcj} . '-' . $seg3;
$loc4=$CHR{$hcj} . '+' . $seg1;
$loc5=$CHR{$hcj} . '+' . $seg2;
$loc6=$CHR{$hcj} . '+' . $seg3;
	foreach $lcj (@{$LOCJ{$loc1}},@{$LOCJ{$loc2}},@{$LOCJ{$loc3}},@{$LOCJ{$loc4}},@{$LOCJ{$loc5}},@{$LOCJ{$loc6}}){   # go over all low conf. junctions
	next unless $L_ATTRIBS{$lcj};	
	if($conf ne "NC"){                          ### added Andre Sep 8th in order to get all exons
		next if $LR{$hcj}*$L_LR{$lcj}>0;    # - require opposite efficiency
	}
	next if $hcj eq $lcj;
	#next unless $STRAND{$hcj} eq $EL_STRAND{$lcj};
	$minDist=min(abs($START{$hcj}-$EL_END{$lcj}),abs($END{$hcj}-$EL_START{$lcj}),abs($START{$hcj}-$EL_START{$lcj}),abs($END{$hcj}-$EL_END{$lcj}) );	
	next unless $minDist<100000;                        # junctions must be closer then 100000 nts
	$lc_jlength=abs($EL_START{$lcj}-$EL_END{$lcj});
	next if $lc_jlength>100000;                         # low conf. junctions must not span more then 100000 nts
	$scale1=($ATTRIBS{$hcj}[1]+$ATTRIBS{$hcj}[2])/2;
	$scale2=($L_ATTRIBS{$lcj}[1]+$L_ATTRIBS{$lcj}[2])/2;
	$LRdist=abs(($ATTRIBS{$hcj}[1]-$ATTRIBS{$hcj}[2])/$scale1+($L_ATTRIBS{$lcj}[1]-$L_ATTRIBS{$lcj}[2])/$scale2 );
	#$RCdist=log ( 1+abs( $ATTRIBS{$hcj}[3]-$ATTRIBS{$hcj}[4] + $L_ATTRIBS{$lcj}[3]-$L_ATTRIBS{$lcj}[4] )  )/log(10);
	$PVdist=abs( log($ATTRIBS{$hcj}[7]+1e-20)/log(10000)-log($L_ATTRIBS{$lcj}[7]+1e-20)/log(10000) );	
	$DSdist=4.5*log($minDist/250+1)/log(100);
	$STRANDdist=6 *  ($STRAND{$hcj} ne $EL_STRAND{$lcj});
	$intsct=intersect( ($START{$hcj},$END{$hcj},$EL_START{$lcj},$EL_END{$lcj}) );
	# goal: from all close enough low confidence jucntions, create a ranking and select the top-one as competing
	# ranking will be done wrt. $Dist (the lower the better)
	# $Dist is a heuristic score which combines several aspects of pairs of junctions, like their minimal distance, are they on the
	# same strand, are bot diff. splice between conditions etc. 
	$Dist=$LRdist+$PVdist+$RCdist+$DSdist+$STRANDdist+ (6*(1-$intsct)*($STRAND{$hcj} eq $EL_STRAND{$lcj}));
	next unless $Dist<$Dist{$hcj};
	$Dist{$hcj}=$Dist;
	$PROXIMAL{$hcj}=$lcj;
	$LRdist{$hcj}=$LRdist;
	$PVdist{$hcj}=$PVdist;
	$RCdist{$hcj}=$RCdist;
	$DSdist{$hcj}=$DSdist;
	$PROX_ABS_DIS{$hcj}=$minDist;
		if ( $minDist==abs($START{$hcj}-$EL_END{$lcj}) ){			
		$PROX_DIS{$hcj}=-$START{$hcj}+$EL_END{$lcj};
		$TYPE{$hcj}="nonME_COMP";		# non-mutual exclusive complex
		}
		if ($minDist==abs($END{$hcj}-$EL_START{$lcj}) ){
		$PROX_DIS{$hcj}=-$END{$hcj}+$EL_START{$lcj};
		$TYPE{$hcj}="nonME_COMP";		
		}

		$TYPE{$hcj}="DUAL_COMP" if ($TYPE{$hcj} eq "nonME_COMP" && $intsct );
		if ($minDist==abs($END{$hcj}-$EL_END{$lcj}) ){
		$PROX_DIS{$hcj}=-$END{$hcj}+$EL_END{$lcj};
		$TYPE{$hcj}="DUAL_COMP";	
		}
		if ($minDist==abs($START{$hcj}-$EL_START{$lcj}) ){
		$PROX_DIS{$hcj}=-$START{$hcj}+$EL_START{$lcj};
		$TYPE{$hcj}="DUAL_COMP";
		}
		
		$CE="no";	#Moved this out of the loop 05-10-14
		if ($END{$hcj}==$EL_END{$lcj} || $START{$hcj}==$EL_START{$lcj}){  # the two junctions have identical start or identical end
		($longj,$shortj)=($hcj,$lcj);
		($longj,$shortj)=($lcj,$hcj) if (abs($END{$hcj}-$START{$hcj})) < (abs($EL_END{$lcj}-$EL_START{$lcj}));
			foreach $ssj (@{$Same_Start{$longj}}) {
			next if abs($EL_START{$longj}-$EL_START{$shortj})<2;	#Added 05-10-14
			next if abs($EL_END{$ssj}-$EL_END{$longj})<2;
			next if ($EL_START{$ssj} eq $EL_END{$shortj}) || ($EL_END{$ssj} eq $EL_START{$shortj});
			next if $ssj eq $shortj || $ssj eq $longj;
			$CE="yes" if ($EL_END{$ssj} < $EL_START{$shortj}-2); 	##  && $L_LR{$ssj} && $L_LR{$ssj}*$L_LR{$shortj}>0);   # both junctions indicate a SE event    # Modified Claudia 05-11-16 to improve annotation
			#print OUT "LOOK1: $hcj\t$lcj\t$longj\t$shortj\t$ssj\t$CE\n";# if $CE eq "yes";
			}
			foreach $sej (@{$Same_End{$longj}}) {
			next if abs($EL_END{$longj}-$EL_END{$shortj})<2;	#Added 05-10-14
			next if abs($EL_START{$sej}-$EL_START{$longj})<2;
			next if ($EL_START{$sej} eq $EL_END{$shortj}) || ($EL_END{$sej} eq $EL_START{$shortj});
			next if $sej eq $shortj || $sej eq $longj;
			$CE="yes" if ($EL_START{$sej} > $EL_END{$shortj}+2); 	##  && $L_LR{$sej}  && $L_LR{$sej}*$L_LR{$shortj}>0)	# Modified Claudia 05-11-16 to improve annotation
			#print OUT "LOOK2: $hcj\t$lcj\t$longj\t$shortj\t$sej\t$CE\n"  if $CE eq "yes";
			}
		}

		# we consider for each high-confidence jucntion several possible competing junctions
		# one pairing might indicate a SE event
		# another one a 3SS or 5SS 
		# we add here all these indications into one competing type
		if ($END{$hcj}==$EL_END{$lcj} ){
		$TYPE{$hcj}="COMP_3'SS_" . ($START{$hcj}-$EL_START{$lcj}) if $STRAND{$hcj} eq "-";
		$TYPE{$hcj}="COMP_5'SS_" . ($START{$hcj}-$EL_START{$lcj}) if $STRAND{$hcj} eq "+";
		$TYPE{$hcj}.="_CE"  if (abs($START{$hcj}-$EL_START{$lcj})>1000 || $CE eq "yes");
		}
		if ($START{$hcj}==$EL_START{$lcj} ){
		$TYPE{$hcj}="COMP_3'SS_" . ($END{$hcj}-$EL_END{$lcj}) if $STRAND{$hcj} eq "+";
		$TYPE{$hcj}="COMP_5'SS_" . ($END{$hcj}-$EL_END{$lcj}) if $STRAND{$hcj} eq "-";
		$TYPE{$hcj}.="_CE" if (abs($END{$hcj}-$EL_END{$lcj})>1000 || $CE eq "yes");
		}
	$TYPE{$hcj}="OPP. STRAND" if $STRAND{$hcj} ne $EL_STRAND{$lcj}
	}
}


if ($conf ne 'IRM'){
print OUT "INCL_COORDs\tGene_Name(s)\tHigh_Confidence_Junction\tCompeting_Junction\tminJDist\tCOMPET_TYPE\tHCJ_5'ss\tHCJ_3'ss\tHCJ_Junc\t";
print OUT "CompJ_5'ss\tCompJ_3'ss\tCompJ_Junc\tHCJ_LR($COND2/$COND1)\tCompJ_LogRatio($COND2/$COND1)\tHCJ_Delta($COND2-$COND1)\tCompJ_Delta($COND2-$COND1)\tHCJ_Pval\t";
print OUT "CompJ_Pval\tHCJ_Eff_$COND1\tHCJ_Eff_$COND2\tHCJ_N_$COND1\tHCJ_N_$COND2\tHCJ_PSI_$COND1\tHCJ_PSI_$COND2\tHCJ_Creads_$COND1\tHCJ_Creads_$COND2\tCompJ_Eff_$COND1\tCompJ_Eff_$COND2\tCompJ_N_$COND1\tCompJ_N_$COND2\tCompJ_PSI_$COND1\tCompJ_PSI_$COND2\tCompJ_Creads_$COND1\tCompJ_Creads_$COND2\n";
}

else {
print OUT "INCL_COORDs\tGene_Name(s)\tHigh_Confidence_Junction\tCOMPET_TYPE\tHCJ_5'ss\tHCJ_3'ss\tHCJ_Junc\t";
print OUT "IRLR\tPvalIR\tHCJ_Delta\tHCJ_Pval\tHCJ_N_$COND1\tHCJ_N_$COND2\tHCJ_PSI_$COND1\tHCJ_PSI_$COND2\n";
}


$,="\t";
foreach $hcj (keys %ATTRIBS){
$lcj=$PROXIMAL{$hcj};
$hcj_gene="";
($NOVD,$NOVA,$NOVJ)=("known","known","known");
($NOVD1,$NOVA1,$NOVJ1)=("NA","NA","NA");
@L_VALS=('NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA');

@mm=split /\_/,$hcj;

$chr=$mm[0];
$strand=$mm[3];
$startj=$chr . '_' . $mm[1] . '_' . $strand;
$min=$mm[1]-1;
$endj=$chr . '_' . $mm[2] . '_'  . $strand;
$max=$mm[2]+1;
$donor=$startj if $strand eq "+";
$donor=$endj if $strand eq "-";
$accpt=$endj if $strand eq "+";
$accpt=$startj if $strand eq "-";
$NOVD="novel" if $DONOR{$donor}<1;
$NOVA="novel" if $ACCPT{$accpt}<1;
$NOVJ="novel" if $JUNCTION{$hcj}<1;

$TYPE{$hcj}.="/RET_INTRON" if  $RET_INTR{$hcj} && ($JUNCTION{$hcj}>0);	#Edited 03-15 from: && (($DONOR{$donor}+$ACCPT{$accpt}<4)||$JUNCTION{$hcj}>0)

unless ($PROXIMAL{$hcj} eq "NA"){
($NOVD1,$NOVA1,$NOVJ1)=("known","known","known");
$PROX_DIS{$hcj}=-$PROX_DIS{$hcj} if $STRAND{$hcj} eq "-";
$lc_jlength=abs($EL_START{$lcj}-$EL_END{$lcj});
$hc_jlength=abs($START{$hcj}-$END{$hcj});
$lcj_gene="";
@l_mm=split /\_/,$PROXIMAL{$hcj};
$min=$l_mm[1]-1 if $l_mm[1]<$mm[1];
$max=$l_mm[2]+1 if $l_mm[2]>$mm[2];
$L_VALS[0]=$L_LR{$lcj};
$L_VALS[1]=$L_DELTA{$lcj};
$L_VALS[2]=$L_ATTRIBS{$lcj}[7];
@L_VALS[3..10]=@{$L_ATTRIBS{$lcj}}[(1..6,8..9)];

@mm1=split /\_/,$PROXIMAL{$hcj};
$chr1=$mm1[0];
$strand1=$mm1[3];
$startj1=$chr1 . '_' . $mm1[1] . '_' . $strand1;
$endj1=$chr1 . '_' . $mm1[2] . '_'  . $strand1;
$donor1=$startj1 if $strand1 eq "+";
$donor1=$endj1 if $strand1 eq "-";
$accpt1=$endj1 if $strand1 eq "+";
$accpt1=$startj1 if $strand1 eq "-";
$NOVD1="novel" if $DONOR{$donor1}<1;
$NOVA1="novel" if $ACCPT{$accpt1}<1;
$NOVJ1="novel" if $JUNCTION{$PROXIMAL{$hcj}}<1;
}





$INCL_COORDS="$chr:$min-$max";
next unless $conf ne 'IRM' || $TYPE{$hcj}=~/RET/;
print OUT "$INCL_COORDS\t";
$,=",";
undef %saw;
@{$Gnames{$hcj}} = grep(!$saw{$_}++, @{$Gnames{$hcj}});
print OUT @{$Gnames{$hcj}};
$,="\t";
if ($conf ne 'IRM'){
print OUT "\t$hcj\t$PROXIMAL{$hcj}\t$PROX_DIS{$hcj}\t$TYPE{$hcj}\t$NOVD\t$NOVA\t$NOVJ\t$NOVD1\t$NOVA1\t$NOVJ1";
print OUT "\t$LR{$hcj}\t$L_VALS[0]\t$DELTA{$hcj}\t$L_VALS[1]\t$ATTRIBS{$hcj}[7]\t$L_VALS[2]\t";
print OUT @{$ATTRIBS{$hcj}}[(1..6,8..9)];
print OUT "\t";
print OUT @L_VALS[3..10];
print OUT "\n";
}

else{
print OUT "\t$hcj\t$TYPE{$hcj}\t$NOVD\t$NOVA\t$NOVJ";
print OUT "\t${$RI{$hcj}}[6]\t${$RI{$hcj}}[7]\t$DELTA{$hcj}\t$ATTRIBS{$hcj}[7]\t";
print OUT @{$ATTRIBS{$hcj}}[3..6];
print OUT "\n";
}

}



close(OUT);



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

sub intersect {	#Return 1 if two segments a,b: @_=(a1,a2,b1,b2) intersect 0 otherwise. No need for @_ to be sorted
my @a = sort { $a <=> $b } @_[0..1];
my @b = sort { $a <=> $b } @_[2..3];
($a_min,$a_max,$b_min,$b_max)=($a[0],$a[1],$b[0],$b[1]);
$isct = !(($a_max < $b_min) || ($b_max < $a_min));
return (0+$isct);
}



