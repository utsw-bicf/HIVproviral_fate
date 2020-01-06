### This is a perl script that,
### Takes the abs.bed file and HIV expression file and gets group number
### Then takes matrix file and changes group number
use strict; use warnings;

my ($HIVexpHigh, $HIVexpMed, $HIVexpLow, $bed, $matrix) = @ARGV;
my (%Highcolumn, %Medcolumn, %Lowcolumn, %Bcolumn, %Ccolumn, %Dcolumn);
my (@Hcol, @medcol, @Lcol, @Macol);

open pdbOUT, ">HIV_added_HighMedLow.pdb" or die $!;

### High
open IN1, "<$HIVexpHigh" or die $!;
while (my $Hlines = <IN1>) {
  chomp $Hlines;
  my @Hcol = split (/\t/, $Hlines);
  push(@{$Highcolumn{$Hcol[0]}},[$Hcol[1]]);
}
close IN1;

open IN2, "<$HIVexpMed" or die $!;
while (my $Hlines = <IN2>) {
  chomp $Hlines;
  my @medcol = split (/\t/, $Hlines);
  push(@{$Medcolumn{$medcol[0]}},[$medcol[1]]);
}
close IN2;

open IN3, "<$HIVexpLow" or die $!;
while (my $Hlines = <IN3>) {
  chomp $Hlines;
  my @Lcol = split (/\t/, $Hlines);
  push(@{$Lowcolumn{$Lcol[0]}},[$Lcol[1]]);
}
close IN3;

open B1, "<$bed" or die $!;
while (my $Blines = <B1>) {
  chomp $Blines;
  my @Bcol = split (/\t/, $Blines);
  if ($Highcolumn{$Bcol[0]}) {
    foreach my $item (@{$Highcolumn{$Bcol[0]}}) {
      if (($item->[0] >= $Bcol[1]) and ($item->[0] <= $Bcol[2])) {
        push(@{$Bcolumn{$Bcol[3]}},[$Bcol[3]]);
      }
    }
  }
  if ($Medcolumn{$Bcol[0]}) {
    foreach my $item (@{$Medcolumn{$Bcol[0]}}) {
      if (($item->[0] >= $Bcol[1]) and ($item->[0] <= $Bcol[2])) {
        push(@{$Ccolumn{$Bcol[3]}},[$Bcol[3]]);
      }
    }
  }
  if ($Lowcolumn{$Bcol[0]}) {
    foreach my $item (@{$Lowcolumn{$Bcol[0]}}) {
      if (($item->[0] >= $Bcol[1]) and ($item->[0] <= $Bcol[2])) {
        push(@{$Dcolumn{$Bcol[3]}},[$Bcol[3]]);
      }
    }
  }
}
close B1;

open M1, "<$matrix" or die $!;
while (my $Mlines = <M1>) {
  chomp $Mlines;
  my @Macol = split(/\s+/, $Mlines);
  if ($Ccolumn{$Macol[1]}) {
    $Mlines =~ s/EDG/ASP/;
  }
  elsif ($Bcolumn{$Macol[1]}) {
    $Mlines =~ s/EDG/ALA/;
  }
  elsif ($Dcolumn{$Macol[1]}) {
    $Mlines =~ s/EDG/SER/;
  }
  print pdbOUT "$Mlines\n";
}
close M1;

