#!/usr/local/bin/perl -w
use strict;
use Getopt::Long;

my (%options, %peaks);
my (@files, @data);
my ($file, $start, $end, $start0, $end0);

GetOptions (\%options, 'dir=s');

opendir(DIR, "$options{'dir'}")||die "Couldn't open the directory: $!\n";
@files = grep { /\.bedpe$/ } readdir(DIR);
close(DIR);

## depending on the strand of mate1 we shift the start/end 4 or 5 bp.
## we take these shifted coordinates as the insertion sites of Tn5
## we report the insertion coordinates as a BED file, which is 0-based for the start, but 1-based for the end

foreach $file (@files){
	print "Shifting $file\n";
	$file =~ /(.+)\.bedpe/;
	open(FILE, "$options{'dir'}/$file")||die "Couldn't open the file $file: $!\n";
	open(OUT, ">$options{'dir'}/$1.insertionSites.bed")||die "Couldn't create the file for $file: $!\n";
	while(<FILE>){
		@data = split(/\t/, $_);
		if($data[8] eq "+"){
			$start = $data[1]+4; $start0 = $start-1;
			$end = $data[5]-5; $end0 = $end-1;
			print OUT "$data[0]\t$start0\t$start\n$data[3]\t$end0\t$end\n";
		}else{
			$end = $data[2]-4; $end0 = $end-1;
			$start = $data[4]+5; $start0 = $start-1;
			print OUT "$data[0]\t$end0\t$end\n$data[3]\t$start0\t$start\n";
		}
	}
	close(FILE);
	close(OUT);
}
