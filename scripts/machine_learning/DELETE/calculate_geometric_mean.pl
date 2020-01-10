### perl script to calculate geometric mean
use strict; use warnings;

##########
my ($counts, $bed) = @ARGV;
my (%column);
my (@col, @col2);
my ($mulct,$ct, $GM);

open IN, "<$counts" or die $!;
while (my $lines = <IN>) {
  chomp $lines;
  my @col = split(/\t/, $lines);
  #push(@{$column{$col[0]}}, [$col[1], $col[2]]);
  #$column->{$col[0]}=$col[1];
  push(@{$column{$col[0]}}, [$col[1]);
}
close IN;

open INF2, "<$bed" or die $!;
while (my $lines2 = <INF2>) {
  chomp $lines2;
  my @col2 = split(/\t/, $lines2);
  $newname = join '_', $col2[0], $col2[1];;
  if ($column{$newname}) {
    $mulct=1;
    $ct=0;
    foreach my $item (@{$column{$col2[0]}}) {
      if (($item->[0] >= $col2[1]) and ($item->[0] <= $col2[1])) {
        $mulct *= $item->[1];
        $ct++;
      } # if
    } # foreach
    $GM = $mulct**(1/$ct);
#    $stat->geometric_mean();
  } # if
  print "$lines2\t$GM\n";
} # while
close INF2;    
