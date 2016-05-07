#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my @args = @ARGV or die "$0 junctionpep_path percolator_path"; # GM18486.rna result_v5/Alu/junctionPep/GM18486.rna.fa result_v5/Alu/cometout/GM18486.rna/0.05

die "No such file found!\n" if scalar(@args) < 1;

(my $sample = basename($args[0])) =~ s/\.fa//g;;
my $fdr = basename($args[1]);

say "# obtain percolator result from '$sample' with fdr '$fdr'";
mkdir 'tmp' if !-e 'tmp';
my @peparr = `awk -F '\t' '\$3 < $fdr' $args[1]/*.pep.percolator | cut -f 5`;
my %pephash;
my %seqhash;
open TMP,">tmp/$sample.tmp";
foreach(0..$#peparr){
	(my $p = $peparr[$_]) =~ s/\R|\.|\*|-//g;
	say TMP ">seq_$_\n",$p if !exists $pephash{$p};
	$pephash{$p} = '';
	$seqhash{"seq_".$_} = $p;
}
close TMP;

say "# formatdb";
system("cp $args[0] tmp/");
system("/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/formatdb -i tmp/$sample.fa -p T");

say "# blast";
my @blastout = `/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/blastall -p blastp -i tmp/$sample.tmp -d tmp/$sample.fa -m8 -F F`;

my %hash = ();
foreach(@blastout){
        chomp;
        my ($id,$uid, $identity,$alignlen) = (split /\s+/)[0,1,2,3];
        my $orilen = length($seqhash{$id});
        if($orilen == $alignlen && $identity >= 100){
		say $seqhash{$id},"\t",$_;
        }
}

system("rm tmp/$sample.*");
