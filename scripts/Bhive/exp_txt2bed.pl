use strict; use warnings;

while(<>) {
  chomp;
    next if ($_ =~ m/^brcd/);
    my ($brcd, $rep, $chr, $coord, $strand, $reads, $mapq, $dna, $rna, $expr) = split(/\t/, $_);
    my ($coord2);
    if ($strand eq '-') {$coord2 = $coord; $coord = $coord2-1;}
    else {$coord2 = $coord+1;}
    print "$chr\t$coord\t$coord2\t$brcd\t$strand\t$rep\t$reads\t$mapq\t$dna\t$rna\t$expr\n";
}
