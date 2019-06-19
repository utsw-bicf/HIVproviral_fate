### This perl program finds the reads in DamIDseq that tough the edges of GATC sites
### perl script
use strict; use warnings;

my ($GATC, $bed) = @ARGV;
my (%column);
my (@col, @col2);

open INF1, "<$GATC" or die $!;
while (my $lines = <INF1>) {
  chomp $lines;
  my @col = split (/\t/, $lines);
  push(@{$column{$col[0]}},[$col[1], $col[2]]);
}

open INF2, "<$bed" or die $!;
while (my $lines2 = <INF2>) {
  chomp $lines2;
  my @col2 = split(/\t/, $lines2);
  if ($column{$col2[0]}) {
    foreach my $item (@{$column{$col2[0]}}) {
      if (($item->[0] - 1) == $col2[2]) {print "$col2[3]\n";}
      elsif (($item->[0] - 1) == $col2[1]) {print "$col2[3]\n";}
      elsif ($item->[1] == $col2[2]) {print "$col2[3]\n";}
      elsif (($item->[1] + 1) == $col2[2]) {print "$col2[3]\n";}
    }
  }
}

