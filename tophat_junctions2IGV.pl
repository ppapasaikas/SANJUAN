
$titleC='track name=junctions description="Diff. Used Junctions CTR"' ."\n";
$titleR='track name=junctions description="Diff. Used Junctions SSA"' ."\n";
open (OUTC,">Diff_Junctions_IGV_Track_CTR.bed");
open (OUTR,">Diff_Junctions_IGV_Track_SSA.bed");
print OUTC $titleC;
print OUTR $titleR;

#open (IN, "temp16") ||die;
$offset=2;	#500 OR 10

while (<>){
$line=$_;
chomp $line;
next unless $line=~/^chr/;


@mat=split /\t/,$line;
@mm1=split /[\_\:]/,$mat[2];
#@mm2=split /[\_\:]/,$mat[2];
$start1=$mm1[1]-$offset;
$end1=$mm1[2]+$offset;
$len1=$end1-$start1-2;
$chr1=$mm1[0];
$strand1=$mm1[3];
$ID1=$mat[2];
$CTRcount1='';
$RNAcount1='';
print OUTC "$chr1\t$start1\t$end1\t$ID1\t$mat[20]\t$strand1\t$start1\t$end1\t255,0,0\t2\t2,2\t0,$len1\n";
print OUTR "$chr1\t$start1\t$end1\t$ID1\t$mat[21]\t$strand1\t$start1\t$end1\t255,0,0\t2\t2,2\t0,$len1\n";


if ($mat[3]=~/^chr/){
@mm2=split /\_/,$mat[3];
$start2=$mm2[1]-$offset;
$end2=$mm2[2]+$offset;
$len2=$end2-$start2-2;
$chr2=$mm2[0];
$strand2=$mm2[3];
$ID2=$mat[3];
$CTRcount2='';
$RNAcount2='';
print OUTC "$chr2\t$start2\t$end2\t$ID2\t$mat[26]\t$strand2\t$start2\t$end2\t0,0,255\t2\t2,2\t0,$len2\n";
print OUTR "$chr2\t$start2\t$end2\t$ID2\t$mat[27]\t$strand2\t$start2\t$end2\t0,0,255\t2\t2,2\t0,$len2\n";
}






}
close IN;


