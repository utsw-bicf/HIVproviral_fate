### perl script
use strict; use warnings;

my ($chrom, $loop) = @ARGV;
my (%column);
my (@col, @col2);
my ($score,$avgscore);
my ($score2, $score3);

open IN, "<$chrom" or die $!;
while (my $lines = <IN>) {
  chomp $lines;
  my @col = split (/\t/, $lines);
  push(@{$column{$col[0]}},[$col[1], $col[2], $col[3]]);
}

open OUT, ">testing_out.txt" or die $!;
open INF2, "<$loop" or die $!;
while (my $lines2 = <INF2>) {
  chomp $lines2;
  my @col2 = split(/\t/, $lines2);
  if ($column{$col2[0]}) {
  $avgscore=0;
  $score3=0;
  $score2=0;
  $score=0;
    foreach my $item (@{$column{$col2[0]}}) {
      if (($item->[0] <= $col2[1]) and ($item->[1] >= $col2[1]) and ($item->[1] <= $col2[2])) {
        $score = $item->[2]*($item->[1] - $col2[1]);
      }
      elsif (($item->[0] >= $col2[1]) and ($item->[0] <= $col2[2]) and ($item->[1] <= $col2[2]) and ($item->[1] >= $col2[1])) {
         $score2 += $item->[2]*($item->[1] - $item->[0]);
      }
      elsif (($item->[0] <= $col2[2]) and ($item->[1] >= $col2[2]) and ($item->[0] >= $col2[1])) {
         $score3 += $item->[2]*($col2[2] - $item->[0]);
      }
      $avgscore = sprintf("%.2f", ($score+$score2+$score3)/($col2[2]-$col2[1]));
    }
print OUT "$lines2\t$avgscore\n";
  }
}
