#!/usr/local/bin/perl -w
use strict;
use Getopt::Long;

my (%options);
my (@files, @stats);
my ($file, $stats, $total);

GetOptions (\%options, 'dir=s');

opendir(DIR, "$options{'dir'}")||die "Couldn't open the directory: $!\n";
@files = grep { /p1\.fq_fastqc\.txt\.gz$/ } readdir(DIR);
close(DIR);

$total = 0;
foreach $file (@files){
	$stats = `zcat $options{'dir'}/$file | grep ^Total`;
	@stats = split(/\s/, $stats);
	$total += $stats[2];
}
$files[0] =~ /(.+)\_unk_somite.+/;
print "$1\t$total\n";
