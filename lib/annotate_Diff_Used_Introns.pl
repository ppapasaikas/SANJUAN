$Diff_Junct_Eff_HC=$ARGV[0];
$Diff_Junct_Eff_LC=$ARGV[1];
$olapSel_Junc2Tx=$ARGV[2];
$Diff_Intr_Ret=$ARGV[3];
$ENSid2Name=$ARGV[4];
$ENS_Tx_Junc=$ARGV[5];
$SuppJun=$ARGV[6]; #Require Supporting Junction switch? N->no
$conf=$ARGV[7];
$COND1=$ARGV[8];
$COND2=$ARGV[9];

($minLFC,$minPvRET)=(0.92,0.0001) if $conf eq 'VHC';	#x3   fold
($minLFC,$minPvRET)=(0.69,0.0010) if $conf eq 'HC';	#x2   fold
($minLFC,$minPvRET)=(0.41,0.0100) if $conf eq 'MC';	#x1.5 fold

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
next unless abs($mat[6]-$mat[5])>$minDPSI;	
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




open (IN, $Diff_Intr_Ret) ||die;
while (<IN>){
$line=$_;
next if $line=~/^>/;
chomp $line;
@mat=split /\t/,$line;
next unless ($mat[5]*$mat[6])<0 || $SuppJun eq 'N';
next if abs($mat[6])<$minLFC;
next unless $mat[7] < $minPvRET;
$RET_INTR{$mat[0]}=defined;
@{$RI{$mat[0]}}=@mat;
}
close IN;


open (IN, $olapSel_Junc2Tx)||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
next unless $RET_INTR{$mat[3]};

	unless ($ENSid2Gname{$mat[9]}){
	${$Gnames{$mat[3]}}[0]="NA";
	next;
	}

	if (${$Gnames{$mat[3]}}[0] eq "NA") {${$Gnames{$mat[3]}}[0]=$ENSid2Gname{$mat[9]}} else {
	push @{$Gnames{$mat[3]}},$ENSid2Gname{$mat[9]}
	}
}
close IN;


print "INCL_COORDs\tGene_Name(s)\tHigh_Confidence_Junction\tCOMPET_TYPE\tHCJ_5'ss\tHCJ_3'ss\tHCJ_Junc\t";
print "IRLR\tPvalIR\tHCJ_Delta\tHCJ_Pval\tHCJ_N_$COND1\tHCJ_N_$COND2\tHCJ_PSI_$COND1\tHCJ_PSI_$COND2\n";


$,="\t";
foreach $hcj (keys %RET_INTR){
$hcj_gene="";
($NOVD,$NOVA,$NOVJ)=("known","known","known");
($NOVD1,$NOVA1,$NOVJ1)=("NA","NA","NA");

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
next unless $JUNCTION{$hcj}>0;	#Edited 03-15 from: && (($DONOR{$donor}+$ACCPT{$accpt}<4)||$JUNCTION{$hcj}>0)
$TYPE{$hcj}.="/RET_INTRON";




$INCL_COORDS="$chr:$min-$max";
print "$INCL_COORDS\t";
$,=",";
undef %saw;
@{$Gnames{$hcj}} = grep(!$saw{$_}++, @{$Gnames{$hcj}});
print @{$Gnames{$hcj}};
$,="\t";
print "\t$hcj\t$TYPE{$hcj}\t$NOVD\t$NOVA\t$NOVJ";
print "\t${$RI{$hcj}}[6]\t${$RI{$hcj}}[7]\t$DELTA{$hcj}\t$ATTRIBS{$hcj}[7]\t";
print @{$ATTRIBS{$hcj}}[3..6];
print "\n";


}
















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



