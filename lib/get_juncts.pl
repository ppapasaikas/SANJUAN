$SuppThreshold=3;	#Minimum Total Number of Junction Reads  
$FractThreshold=0.6;	#Minimum Fraction of Replicates with Junction
$SAMFILE=$ARGV[0];
$low_seq_req=$ARGV[1]; # Y or N
$fn_out=$ARGV[2];

open (OUT,">".$fn_out) || die $!;

open (IN, "$SAMFILE") ||die;	#Open Merged Junctions File

while (<IN>){
$line=$_;
@mat=split /\t/,$line;


#next unless length($mat[2])<8;	#Remove junctions coming from non-canonical chromosomes. 03-2015
next if ( (length($mat[2])>7 && $mat[2]=~/chr/) ||  (length($mat[2])>4 && $mat[2]!~/chr/) );

next unless $mat[5]=~/^[\dMN]+$/;

@CIG=split /[MN]/,$mat[5];	#Split CIGAR field in spanned segments (e.g 30M400N50M ->30, 400, 50)
$nj=int($#CIG/2);		#For each two additional fields (#N#M) add one junction.

$strand="NA";
#83,99,147 and 163.
if ($mat[1]==99 || $mat[1]==147 || $mat[1]==97 || $mat[1]==145 ){$strand="-"}		#only correct for firstrand
elsif ($mat[1]==83 || $mat[1]==163 || $mat[1]==81 || $mat[1]==161){$strand="+"}		#Only correct for firstrand
if($low_seq_req eq "N" && $strand eq "NA"){next;}
$strand="-" if $line=~/XS:A:\-/;		#Added on 03-2015
$strand="+" if $line=~/XS:A:\+/;		#Added on 03-2015
@RID=split /:/,$mat[0];
$nfields=$#RID-3;
$flowcell=join ('_',@RID[0..$nfields]);
$flowcells++ unless $SAW{$flowcell};
$SAW{$flowcell}=1;
	for $n (1..$nj){
	$st=$mat[3] + $CIG[0] + ($CIG[1] + $CIG[2])*($n>1) + ($CIG[3] + $CIG[4])*($n>2) -1;	#Junction Read Start
	$en=$mat[3] + $CIG[0] + $CIG[1]  + ($CIG[2] + $CIG[3])*($n>1) + ($CIG[4] + $CIG[5])*($n>2) -1;	#Junction Read End
	next if $CIG[0+2*($n-1)]<2;	#Minimum overhang length. Need to correct, does not take into account D\I
	next if $CIG[2+2*($n-1)]<2;
	$id=$mat[2] . '_' . $st  . '_' . $en  . '_' . $strand;
	$fcount{$id}++ unless $JSAW{$id}{$flowcell};
	$JSAW{$id}{$flowcell}=1;
	$count{$id}++;
	}
}


foreach $id (keys %count){
@INF=split /_/,$id;
$opps=$INF[3];
$opps=~ tr /-\+/\+-/;#\
$oppid=$INF[0] . '_' . $INF[1]  . '_' . $INF[2]  . '_' . $opps;
next if $count{$oppid}>$count{$id};	#Spurious junction on opposite strand
if($low_seq_req eq "N" && $fcount{$id}< $FractThreshold * $flowcells){next;} #Absent in significant fraction of flowcells
next if $count{$id}< $SuppThreshold;	#Low count
print OUT "$INF[0]\t$INF[1]\t$INF[2]\t$id\t$count{$id}\t$INF[3]\n";#$fcount{$id}\t$flowcells\n";
}


close(IN);
close(OUT);

# wait 5 min so that output on cluster is written completely for sure 
#sleep 60;









