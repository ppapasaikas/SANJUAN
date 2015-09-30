#open (IN,"/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/UCSC_ENSEMBL_hg19/UCSC_ENSEMBL_hg19.gtf")||die;
#open (IN,"/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/GENCODE19/gencode.v19.annotation.gtf")||die;
open (IN,"/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/danRer10_ENSEMBL/danRer10_GRCz10.81.gtf")||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
next unless $mat[2]=~/exon/i;
$TxID=$1 if $line=~/transcript_id \"([^\"]+)\"/;
$Gene{$TxID}=$1 if $line=~/gene_name \"([^\"]+)\"/;
push @{$STARTS{$TxID}},$mat[3];
push @{$ENDS{$TxID}},$mat[4];
$CHR{$TxID}=$mat[0];
$STRAND{$TxID}=$mat[6];
}
close IN;

foreach $TxID (keys %STARTS) {
$TxStart{$TxID}=min(@{$STARTS{$TxID}});
$TxEnd{$TxID}=max(@{$ENDS{$TxID}});

#print "$CHR{$TxID}\t$TxStart{$TxID}\t$TxEnd{$TxID}\t$TxID\t0\t$STRAND{$TxID}\n";	#Generate Bedfile
###print "$Gene{$TxID}\n";
}
#die;


foreach $TxID (keys %STARTS) {
@S = sort { $a <=> $b } @{$STARTS{$TxID}};
@E = sort { $a <=> $b } @{$ENDS{$TxID}};
	for $i (0..$#S-1){
	$end=$S[$i+1]-1;	#To match the used junction convention
	print $TxID ."\t" . $CHR{$TxID} . '_' . $E[$i]  . '_' . $end  . '_' . $STRAND{$TxID} . "\n";	#Generate Tx2Jn table
	}

}
die;








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
