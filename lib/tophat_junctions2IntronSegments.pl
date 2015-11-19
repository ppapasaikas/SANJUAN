#Use this script to generate k bed files with each one coveringh 1/k length of Intronic Segment:

$JUNCTIONS=$ARGV[0];	#ALL2_junctions.bed

open (IN, $JUNCTIONS) ||die;
while (<IN>){
$line=$_;
chomp $line;
next if $line=~/^tr/;
@mat=split /\t/,$line;
@mm=split /\,/,$mat[10];
@mm2=split /\,/,$mat[11];

$start=$mat[1]+$mm[0];
$end=$mat[2]-$mm[1];

$length=abs($end-$start);
$k=3+int (0.5+$length/1000);	 #Previously k set to 5. Now Dynamic in the range 2->8 (0-150 ->2, 150-500 -> 3, 500-1500->4, 1500-2500->5, 2500-3500->5, >3500 -> 7)
$k=7 if $k>7;
$k=2 if $length<150;

$ktile=int (abs( ($end-$start)/$k ) );

$id=$mat[0] . "_" . $start . "_" . $end . "_" . $mat[5];
next if $saw{$id};
$saw{$id}=defined;

$COUNT{$id}++;
$SUPP{$id}+=$mat[4];
$STRAND{$id}=$mat[5];

for (1..$k){
$st=$start+($_-1)*$ktile;
$en=$start+($_)*$ktile;
$st=$start+1 if $_==1;
$en=$end-1 if $_==$k;
next if ($en-$st)<5;
$sid=$id . "_" . $_;
$BED="$mat[0]\t$st\t$en\t$sid\t0\t$STRAND{$id}";
print $BED . "\n";
}
}






