#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;

my (%options);
my ($dir, $total, $unique, $multi, $unmapped);

GetOptions (\%options, 'dir=s');

my $path = $options{'dir'};

print "sample\ttotal\tunique\tmultimapped\tunmapped(%)\n";
opendir(DIR, $path)||die "Could not open the directory $options{'dir'}: $!\n";
while($dir = readdir(DIR)){ 
    $total = 0;    $unique = 0;    $multi = 0;    $unmapped = 0;
    open(FILE, "$path\/$dir\/Log.final.out")||next;
    while(<FILE>){
        if($_ =~ /Number of input reads \|\t(\d+)$/){ $total = $1; }
        if($_ =~ /Uniquely mapped reads number \|\t(\d+)$/){ $unique = $1; }
        if($_ =~ /Number of reads mapped to multiple loci \|\t(\d+)$/){ $multi = $1; }
	      if($_ =~ /Number of reads mapped to too many loci \|\t(\d+)$/){ $multi += $1; }
        if($_ =~ /.+unmapped.+ \|\t(\d+\.\d+)\%/){ $unmapped += $1; }
    }
    close(FILE);
    print "$dir\t$total\t$unique\t\t$multi\t\t$unmapped\n";
}
closedir(DIR);
