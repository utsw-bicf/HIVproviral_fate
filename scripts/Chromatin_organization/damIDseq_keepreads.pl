### This perl program finds the reads in DamIDseq that touch the edges of GATC sites
### perl script
use strict; use warnings;

my ($GATC, $bed) = @ARGV;
my $column={};
my $column2={};
my (@col, @col2);

open INF1, "<$GATC" or die $!;
while (my $lines = <INF1>) {
  chomp $lines;
  my @col = split (/\t/, $lines);
  $column->{$col[0]}->{$col[1]-1}=1;
  $column2->{$col[0]}->{$col[2]}=1;
}

open INF2, "<$bed" or die $!;
while (my $lines2 = <INF2>) {
  chomp $lines2;
  my @col2 = split(/\t/, $lines2);
  if ($column->{$col2[0]}->{$col2[2]}) {print "$col2[3]\n";}
  elsif ($column2->{$col2[0]}->{$col2[1]}) {print "$col2[3]\n";}
}
