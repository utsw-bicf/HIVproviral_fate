use strict; use warnings;

while(<>) {
  chomp;
    next if ($_ =~ m/^brcd/);
    my ($brcd, $chr, $coord, $strand, $reads, $mapq, $rep) = split(/\t/, $_);
    my ($coord2);
    if ($strand eq '-') {$coord2 = $coord; $coord = $coord2-1;}
    else {$coord2 = $coord+1;}
    print "$chr\t$coord\t$coord2\t$brcd\t$strand\t$reads\t$mapq\t$rep\n";
}
