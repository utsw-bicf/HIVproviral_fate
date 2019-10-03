### perl script
use strict; use warnings;

my ($damID, $hivbed) = @ARGV;
my (%column);
my (@col, @col2);
my ($score,$score2, $score3, $score4, $totscore);

open IN, "<$damID" or die $!;
while (my $lines = <IN>) {
  chomp $lines;
  my @col = split (/\t/, $lines);
  push(@{$column{$col[0]}},[$col[1], $col[2], $col[3]]);
}
close IN;

open INF2, "<$hivbed" or die $!;
while (my $lines2 = <INF2>) {
  chomp $lines2;
  my @col2 = split(/\t/, $lines2);
  if ($column{$col2[0]}) {
  $totscore=0;
  $score3=0;
  $score2=0;
  $score=0;
  $score4=0;
    foreach my $item (@{$column{$col2[0]}}) {
      if (($item->[0] <= $col2[1]) and ($item->[1] >= $col2[1]) and ($item->[1] <= $col2[2])) {
        $score = $item->[2]/((1 + $item->[1] - $col2[1])/($item->[1] - $item->[0] + 1));
      }
      elsif (($item->[0] >= $col2[1]) and ($item->[0] <= $col2[2]) and ($item->[1] <= $col2[2]) and ($item->[1] >= $col2[1])) {
        $score2 += $item->[2];
      }
      elsif (($item->[0] <= $col2[2]) and ($item->[1] >= $col2[2]) and ($item->[0] >= $col2[1])) {
        $score3 += $item->[2]/(($col2[2]- $item->[0] + 1)/($item->[1] - $item->[0] + 1));
      }
      elsif (($item->[0] <= $col2[1]) and ($item->[1] >= $col2[2])) {
        $score4 = $item->[2]/((($col2[2] - $col2[1]) + 1)/($item->[1] - $item->[0] + 1));
      }
      $totscore = sprintf("%.2f", ($score + $score2 + $score3 + $score4));
    }
print "$lines2\t$totscore\n";
  }
}
close INF2;
