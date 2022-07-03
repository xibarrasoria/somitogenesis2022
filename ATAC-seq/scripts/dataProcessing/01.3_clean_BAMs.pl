#!/usr/local/bin/perl -w
use strict;
use Getopt::Long;

my (%options);
my (@samples);
my ($dir, $sample);

GetOptions (\%options, 'file=s');

## Read params
open(FILE, "$options{'file'}")||die "Couldn't open the input file: $!\n";
while(<FILE>){
	if($_ =~ /^dir\s*=\s*(.+)/){ $dir = $1; }
	if($_ =~ /^samples\s*=\s*(.+)/){ @samples = split(/\,/, $1); }
}
close(FILE);

## Prepare the output directory
system("mkdir $dir/DUPSmetrics");

## Write out commands
open(OUT, ">01.3_clean_BAMs.sh")||die " Couldn't create output file: $!\n";
foreach $sample(@samples){
	## remove PCR duplicates
	print OUT "java -jar \$picard MarkDuplicates I=$dir/BWA/$sample.bam O=$dir/BWA/$sample.noDUPs.bam M=$dir/DUPSmetrics/$sample.markDup_metrics.txt REMOVE_DUPLICATES=TRUE\n";

	## remove low quality alignments, supplementary alignments, singletons and alignments not in the autosomes or the X chr
	print OUT "samtools view -bh -q 30 -f 0x02 -F 0x800 -L $dir/mm10_usedChr.bed $dir/BWA/$sample.noDUPs.bam > $dir/BWA/$sample.noDUPs.GQ.bam\n";
	
	## index
	print OUT "samtools index $dir/BWA/$sample.noDUPs.GQ.bam\n";
	
	## final library size
	print OUT "samtools flagstat $dir/BWA/$sample.noDUPs.GQ.bam >> $dir/BWA/$sample.flagstat\n";

	## remove intermediate files
	print OUT "rm $dir/BWA/$sample.bam\n";
	print OUT "rm $dir/BWA/$sample.noDUPs.bam\n\n";
}
close(OUT);
