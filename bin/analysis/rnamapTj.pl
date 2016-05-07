#! /u/home/y/ybwang/perl
use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my @dir = glob '/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/RNA/*';

foreach(@dir){
	my $log = $_.'/Log.final.out';
	my @r = `grep 'Uniquely mapped reads number' $log | cut -f 2`;
	(my $f = basename($_)) =~ s/\.rna//g;
	print $f,"\t",$r[0];
	
}
