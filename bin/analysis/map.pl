#! /u/home/y/ybwang/perl
use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my $a = 'map.txt';
my $b = 'result_v5/Alu.0.05.list';

my %hash;
open IN, $a;
while(<IN>){
	chomp;
	my @arr = split /\t/;
	$hash{$arr[3]} = '';
}
close IN;

open AS,$b;
while(<AS>){
        chomp;
	if(!/^chr/){
		say;
	}else{
		say if exists $hash{$_};
	}
}
close AS;
