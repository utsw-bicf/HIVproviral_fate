use strict;
use warnings;

open FILE1, "/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/SRR5261760_GATC_counts.bed";
open FILE2, "/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/SRR5261762_GATC_counts.bed";
open INPUT1, "/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/SRR5261759_GATC_counts.bed";
open INPUT2, "/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/SRR5261761_GATC_counts.bed";
open FILECAT, "/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/cat_beds_sorted.bed";

my %count_table;
my $ID;

### Added this next while, which is the concatenation of all files to get the names first
while (my $line = <FILECAT>) {
	chomp $line;
        next if ($line !~ m/^chr[1-9YXM]/);
	my @tmp = split(/\t/, $line);
	my $chr = $tmp[0];
	my $start = $tmp[1];
	$ID = join '_', $tmp[3], $tmp[0], $tmp[1];
}


while (my $line = <FILE1>) {
	chomp $line;
next if ($line !~ m/^chr[1-9YXM]/);
	my @tmp = split(/\t/, $line);
	my $chr = $tmp[0];
	my $start = $tmp[1];
	$ID = join '_', $tmp[3], $tmp[0], $tmp[1];
	my $count = $tmp[6];
	$count_table{$ID}=$count;
}

while (my $line = <FILE2>) {
	chomp $line;
next if ($line !~ m/^chr[1-9YXM]/);
	my @tmp = split(/\t/, $line);
	my $chr = $tmp[0];
	my $start = $tmp[1];
	$ID = join '_', $tmp[3], $tmp[0], $tmp[1];
	my $count = $tmp[6];
	$count_table{$ID}=$count_table{$ID}."\t".$count;
}

while (my $line = <INPUT1>) {
	chomp $line;
next if ($line !~ m/^chr[1-9YXM]/);
	my @tmp = split(/\t/, $line);
	my $chr = $tmp[0];
	my $start = $tmp[1];
	$ID = join '_', $tmp[3], $tmp[0], $tmp[1];
	my $count = $tmp[6];
	$count_table{$ID}=$count_table{$ID}."\t".$count;
}

while (my $line = <INPUT2>) {
	chomp $line;
next if ($line !~ m/^chr[1-9YXM]/);
	my @tmp = split(/\t/, $line);
	my $chr = $tmp[0];
	my $start = $tmp[1];
	$ID = join '_', $tmp[3], $tmp[0], $tmp[1];
	my $count = $tmp[6];
	$count_table{$ID}=$count_table{$ID}."\t".$count;
}


close FILE1;
close FILE2;
close INPUT1;
close INPUT2;


open NEWFILE, ">/project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/damIDseq/damIDseq_sbs/DE/test.tsv";

print NEWFILE "frag_id\tLaminB11\tLaminB12\tDam1\tDam2\n";
foreach my $ID (keys %count_table) {
        next if ($ID !~ m/^chr[1-9YXM]/);
	print NEWFILE "$ID\t$count_table{$ID}\n";
}

close NEWFILE;
