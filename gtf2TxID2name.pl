open (IN,"/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/danRer10_ENSEMBL/danRer10_GRCz10.81.gtf")||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
$GeneName="";
next unless $mat[2]=~/transcript/i;
next unless $mat[8]=~/transcript_id \"([a-zA-Z0-9]+)\";/;
$TxID=$1;
next unless $mat[8]=~/gene_id \"([a-zA-Z0-9]+)\";/;
$GeneID=$1;
$GeneName=$1 if $mat[8]=~/gene_name \"([^\"]+)\";/;
$GeneName=$GeneID if length($GeneName)<2;
print "$TxID\t$GeneName\n";
}
