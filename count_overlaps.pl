#$report_events="N";	#Should the overlapping events be printed out? (Y/N)
$exclude="N";		#Is there a third file containing events that must be excluded? (Y/N)

$File1=$ARGV[0];
$File2=$ARGV[1];
$OUTFile=$ARGV[2];
$ExclusionFile=$ARGV[3];

$helpful_message="\nUsage: perl count_overlaps.pl File1 File2 OutFile ControlFile\n
File1, File2: SANJUAN output files to be compared.\n
OutFile: Output file where the overlapping events will be written.\n
ControlFile: Optional SANJUAN output file containing events to be excluded from the comparison.\n";

if ($#ARGV<2 || $#ARGV>3){
die "$helpful_message";
}
$exclude="Y" if $#ARGV==3;



#####################################################################
#####################################################################


if ($exclude eq "Y"){
open (EXCL, $ExclusionFile)||die "$!";
while (<EXCL>){
$line=$_;
next if $line=~/^INCL/;
@mat=split /\t/,$line;
$ID=$mat[0];
$exclude{$ID}=defined;
}
}




open (IN1, $File1)||die "$!";

while (<IN1>){
$line=$_;
$title=$line if $line=~/^INCL/;
next if $line=~/^INCL/;
@mat=split /\t/,$line;
$ID=$mat[0];
$INFO{$ID}=$line;
next if $exclude{$ID} && $exclude eq "Y";
next if $mat[5] eq "UNIDNT_COMP";
next if $mat[5] eq "nonME_COMP";
next if ($event1{$ID});
$event1{$ID}=defined;
$c1++;
$c1_3SS++ if $mat[5]=~/^COMP_3'SS_[\-\d]+$/;
$c1_5SS++ if $mat[5]=~/^COMP_5'SS_[\-\d]+$/;
$c1_CE++ if $mat[5]=~/^COMP_[35]'SS_[\-\d]+_CE$/;
$c1_RI++ if $mat[5] eq 'UNIDNT_COMP/RET_INTRON';
$TYPE{$ID}=$mat[5];
$JUNC{$ID}=$mat[8].$mat[11];
$SPLS{$ID}=$mat[6].$mat[7].$mat[9].$mat[10];
next if $SPLS{$ID}=~/NA/;
$C1_NOVEL++ if $JUNC{$ID}=~/novel/;
$C1_KNOWN++ if $JUNC{$ID}!~/novel/;
$C1_CRYPT++ if $SPLS{$ID}=~/novel/;
}

open (IN2, $File2)|| die "$!";
open (OUT, ">$OUTFile") || die "$!";

while (<IN2>){
$line=$_;
next if $line=~/^INCL/;
@mat=split /\t/,$line;

$ID=$mat[0];
next if $exclude{$ID} && $exclude eq "Y";
next if $mat[5] eq "UNIDNT_COMP";
next if $mat[5] eq "nonME_COMP";
next if ($event2{$ID});
$event2{$ID}=defined;
$c2++;
$c2_3SS++ if $mat[5]=~/^COMP_3'SS_[\-\d]+$/;
$c2_5SS++ if $mat[5]=~/^COMP_5'SS_[\-\d]+$/;
$c2_CE++ if $mat[5]=~/^COMP_[35]'SS_[\-\d]+_CE$/;
$c2_RI++ if $mat[5] eq 'UNIDNT_COMP/RET_INTRON';
$TYPE{$ID}=$mat[5];
$JUNC{$ID}=$mat[8].$mat[11];
$SPLS{$ID}=$mat[6].$mat[7].$mat[9].$mat[10];
next if $SPLS{$ID}=~/NA/;
$C2_NOVEL++ if $JUNC{$ID}=~/novel/;
$C2_KNOWN++ if $JUNC{$ID}!~/novel/;
$C2_CRYPT++ if $SPLS{$ID}=~/novel/;
}

print OUT $title;
foreach (keys %event1){
next unless $event2{$_};
$c++;
print OUT $INFO{$_};
$c_3SS++ if $TYPE{$_}=~/^COMP_3'SS_[\-\d]+$/;
$c_5SS++ if $TYPE{$_}=~/^COMP_5'SS_[\-\d]+$/;
$c_CE++ if  $TYPE{$_}=~/^COMP_[35]'SS_[\-\d]+_CE$/;
$c_RI++ if  $TYPE{$_} eq 'UNIDNT_COMP/RET_INTRON';
next if $SPLS{$_}=~/NA/;
$C_NOVEL++ if $JUNC{$_}=~/novel/;
$C_KNOWN++ if $JUNC{$_}!~/novel/;
$C_CRYPT++ if $SPLS{$_}=~/novel/;
}
close IN2;
close OUT;

#die "Done!\n" if $report_events eq "Y";

print "TOTAL_SIMPLE_EVENTS:\t$c1\t$c2\t$c\n";
print "SIMPLE KNOWN EVENTS:\t$C1_KNOWN\t$C2_KNOWN\t$C_KNOWN\n";
print "SIMPLE NOVEL EVENTS:\t$C1_NOVEL\t$C2_NOVEL\t$C_NOVEL\n";
print "SIMPLE CRYPTIC ss:\t$C1_CRYPT\t$C2_CRYPT\t$C_CRYPT\n";
print "TOTAL_3'SS_EVENTS:\t$c1_3SS\t$c2_3SS\t$c_3SS\n";
print "TOTAL_5'SS_EVENTS:\t$c1_5SS\t$c2_5SS\t$c_5SS\n";
print "TOTAL_CE_EVENTS:\t$c1_CE\t$c2_CE\t$c_CE\n";
print "TOTAL_RI_EVENTS:\t$c1_RI\t$c2_RI\t$c_RI\n";


