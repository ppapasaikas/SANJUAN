while (<>){
$line=$_;
if ($line=~/^INCL/){
print $line;
next;
}
chomp $line;
@mat=split /\t/,$line;
$ID1=$mat[2] . '_' .  $mat[3];
$IDrev{$ID1}=$mat[3] . '_' .  $mat[2];
next if $mat[5] eq "UNIDNT_COMP";
next if $mat[5] eq "nonME_COMP";
$SAW_EVENT{$ID1}=1;

@{$INFO{$ID1}}=@mat;

}


$,="\t";
foreach $ID1 (keys %INFO){
($L1,$L2)=(0,0);
$ID=$ID1;
$ID2=$IDrev{$ID1};
next if ($WROTE_EVENT{$ID1});
next if ($WROTE_EVENT{$ID2});
	if ($SAW_EVENT{$ID2}){
	@m1=split /_/,$ID1;
	@m2=split /_/,$ID2;
	$L1=abs($m1[2]-$m1[1]);
	$L2=abs($m2[2]-$m2[1]);
	$ID=$ID2 if $L2<$L1;
	}
$WROTE_EVENT{$ID1}=1;
$WROTE_EVENT{$ID2}=1;

print @{$INFO{$ID}};
print "\n";
}

