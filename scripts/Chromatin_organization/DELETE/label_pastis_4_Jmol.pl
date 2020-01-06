### This is a perl script that,
### Takes the abs.bed file and HIV expression file and gets group number
### Then takes matrix file and changes group number
use strict; use warnings;

my ($HIVexp, $bed, $matrix) = @ARGV;
my (%Hcolumn, %Bcolumn);
my (@Hcol, @Mcol);

open OUT, ">testing_out.txt" or die $!;
open pdbOUT, ">HIV_added.pdb" or die $!;

open IN, "<$HIVexp" or die $!;
while (my $Hlines = <IN>) {
  chomp $Hlines;
  my @Hcol = split (/\t/, $Hlines);
  push(@{$Hcolumn{$Hcol[0]}},[$Hcol[1]]);
}
close IN;

open INF2, "<$bed" or die $!;
while (my $Blines = <INF2>) {
  chomp $Blines;
  my @Bcol = split (/\t/, $Blines);
  if ($Hcolumn{$Bcol[0]}) {
    foreach my $item (@{$Hcolumn{$Bcol[0]}}) {
      if (($item->[0] >= $Bcol[1]) and ($item->[0] <= $Bcol[2])) {
        push(@{$Bcolumn{$Bcol[3]}},[$Bcol[3]]);
        print OUT "$Bcol[3]\n";
      }
    }
  }
}
close INF2;

open INF3, "<$matrix" or die $!;
while (my $Mlines = <INF3>) {
  chomp $Mlines;
  my @Mcol = split(/\s+/, $Mlines);
  if ($Bcolumn{$Mcol[1]}) {
    $Mlines =~ s/EDG/ALA/;
  }
  print pdbOUT "$Mlines\n";
}
