#! /u/home/y/ybwang/perl
use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my @dir = glob '/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/RNA/*';

my $cometdir = '/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/result_v5/Alu/cometout/PJ_log_out/';

foreach(@dir){
	(my $f = basename($_)) =~ s/\.rna//g;
	my @grep = ` grep -c $f $cometdir/*`;
	my @arr = ();
	foreach my $r (@grep){
		chomp($r);	
		my $n = (split /:/,$r)[1];
		push @arr, basename((split /:/,$r)[0]) if $n > 0;
	}
	next if scalar(@arr) < 1;
	my $sum = 0;
	foreach $a (@arr){
		my @loadgrep = `grep 'Load spectra:' $cometdir/$a | cut -d ':' -f 2 | sed 's/ //g' | awk 'BEGIN{sum=0}{sum+=\$_}END{print sum}'`;
		my $n = $loadgrep[0];
		chomp($n);
		$sum += $n;
	}
	say $f,"\t",$sum;
}
