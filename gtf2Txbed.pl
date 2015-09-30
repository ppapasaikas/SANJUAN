open (IN,"/users/jvalcarcel/ppapasaikas/TOPHAT_INDEXES/danRer10_ENSEMBL/danRer10_GRCz10.81.gtf")||die;
while (<IN>){
$line=$_;
chomp $line;
@mat=split /\t/,$line;
next unless $mat[2]=~/transcript/i;
next unless $mat[8]=~/transcript_id \"([a-zA-Z0-9]+)\";/;
$TxID=$1;
print "$mat[0]\t$mat[3]\t$mat[4]\t$TxID\t0\t$mat[6]\n";
}
