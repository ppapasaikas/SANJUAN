open (IN1,"sanjuan_KD1/Annotated_Diff_Junctions.txt")||die;
while (<IN1>){
$line=$_;
print $line if $line=~/^INCL/;
chomp $line;
next if $line=~/^INCL/;
@mat=split /\t/,$line;
$ID=$mat[0];

next if $mat[5] eq "UNIDNT_COMP";
next if $mat[5] eq "nonME_COMP";
next if ($event1{$ID});
$event1{$ID}=defined;
@{$INFO1{$ID}}=@mat;
$c1++;
$c1_3SS++ if $mat[5]=~/^COMP_3'SS_[\-\d]+$/;
$c1_5SS++ if $mat[5]=~/^COMP_5'SS_[\-\d]+$/;
$c1_CE++ if $mat[5]=~/^COMP_[35]'SS_[\-\d]+_CE$/;
$c1_RI++ if $mat[5] eq 'UNIDNT_COMP/RET_INTRON';
}


open (IN2,"sanjuan_KD2/Annotated_Diff_Junctions.txt")||die;
while (<IN2>){
$line=$_;
chomp $line;
next if $line=~/^INCL/;
@mat=split /\t/,$line;
$ID=$mat[0];
next if $mat[5] eq "UNIDNT_COMP";
next if $mat[5] eq "nonME_COMP";
next if ($event2{$ID});
$event2{$ID}=defined;
@{$INFO2{$ID}}=@mat;
$c2++;
$c2_3SS++ if $mat[5]=~/^COMP_3'SS_[\-\d]+$/;
$c2_5SS++ if $mat[5]=~/^COMP_5'SS_[\-\d]+$/;
$c2_CE++ if $mat[5]=~/^COMP_[35]'SS_[\-\d]+_CE$/;
$c2_RI++ if $mat[5] eq 'UNIDNT_COMP/RET_INTRON';
$TYPE{$ID}=$mat[5];
}


$,="\t";
foreach (keys %event1){
next unless $event2{$_};
print @{$INFO1{$_}};
print "\n";
}



