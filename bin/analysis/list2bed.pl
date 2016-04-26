#! /u/home/y/ybwang/perl
use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my $in = '1.region'; #'result_v5/Alu.0.05.list';

my $head = '';
open IN, $in;
while(<IN>){
	chomp;
	if(/======/){
		($head = $_) =~ s/======|\.rna//g;
	}
	next if !/chr/;
	my @arr = split /,|_/;
	my $id = $arr[0]."_".$arr[1]."_".$arr[2]."_".$arr[11]."_".$head;
	say $arr[3],"\t",$arr[5]-1,"\t",$arr[7],"\t",$_,"\t0\t",$arr[4];
	say $arr[3],"\t",$arr[8]-1,"\t",$arr[10],"\t",$_,"\t0\t",$arr[4];
	
}
