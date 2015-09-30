#open (IN, "temp33")||die;

while (<>){
$line=$_;	
@mat=split /\t+/,$line;

if ($mat[5]=~/3'SS/ && $mat[5]!~/CE/) {
$countA3SS++ unless $A3SS{$mat[2]} || $A3SS{$mat[3]};
$A3SS{$mat[2]}=defined;
$A3SS{$mat[3]}=defined;
}

if ($mat[5]=~/5'SS/ && $mat[5]!~/CE/) {
$countA5SS++ unless $A5SS{$mat[2]} || $A5SS{$mat[3]};
$A5SS{$mat[2]}=defined;
$A5SS{$mat[3]}=defined;
}

if ($mat[5]=~/_CE/) {
$countCE++ unless $CE{$mat[2]} || $CE{$mat[3]};
$CE{$mat[2]}=defined;
$CE{$mat[3]}=defined;
}	
	
if ($mat[5]=~/RET/) {
$countRI++ unless $RI{$mat[2]};
}
	
}
print "ALT3'ss\tALT5'ss\tCEx\tRIs\n";
print "$countA3SS\t$countA5SS\t$countCE\t$countRI\n";
