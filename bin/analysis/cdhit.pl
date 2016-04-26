#! /u/home/y/ybwang/perl
use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my $in = 'result_v5/Alu.0.05.list';

# read fasta file
my $seq = readSeq();

# read list file and write each junctino seq
my %junction;
my $head = '';
open IN, $in;
while(<IN>){
	chomp;
	if(/======/){
		($head = $_) =~ s/======|\.rna//g;
	}
	next if !/chr/;
	my @arr = split /,|_/;
	my $id = join "_", @arr[0,1,2,4];
	$junction{$id}{$_."_".$head} = $seq->{$_."_".$head};
}
close IN;

mkdir 'tmp' if !-e 'tmp';
foreach my $k1(keys %junction){
	next if scalar(keys %{$junction{$k1}}) < 2; # output junction above 2
	open OUT, ">tmp/".$k1;
	foreach my $k2(keys %{$junction{$k1}}){
		say OUT ">",$k2,"\n",$junction{$k1}{$k2};
	}
	close OUT;
	# cdhit search
	cdhit("tmp/".$k1);
}

# --sub--
sub readSeq {
	my @list = glob "result_v5/Alu/junctionPep/*.fa";
	my %seq = ();
	foreach(@list){
		(my $basename = basename($_)) =~ s/\..*$//g;
		open IN, $_;
		my $id = '';
		while(<IN>){
			chomp;
			if(/>/){
				($id = $_) =~ s/>//g;
			}else{
				$seq{$id."_".$basename} = $_;
			}
		}
		close IN;
	}
	return \%seq;
}

sub cdhit {
	my $in = shift;
	my $out = $in.'.cdhit';
	my $bin = '/u/home/y/ybwang/nobackup-yxing/program/cd-hit-v4.6.5-2016-0304/cd-hit';
	say "$bin -i $in -d 0 -o $out -c 0.9 -n 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0";
	system("$bin -i $in -d 0 -o $out -c 0.9 -n 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0");	
}
