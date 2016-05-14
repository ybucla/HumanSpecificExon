#! /u/home/y/ybwang/perl
# merge junction sequence for single individuals with cdhit (used for human tissue proteome)
use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my $in = shift; # 'liver'
my $fasta = "single_junctionPep/$in*";

# read fasta file
my $seq = readSeq($fasta);

# read list file and write each junctino seq
my %junction;
my $head = '';
foreach my $head(keys %{$seq}){
	my @arr = split /,|_/,$head;
        my $id = join "_", @arr[0,1,2,4];
	say $id;
        $junction{$id}{$head} = $seq->{$head};
}

mkdir 'tmp' if !-e 'tmp';
foreach my $k1(keys %junction){
	#next if scalar(keys %{$junction{$k1}}) < 2; # output junction above 2
	open OUT, ">tmp/".$k1;
	foreach my $k2(keys %{$junction{$k1}}){
		say OUT ">",$k2,"\n",$junction{$k1}{$k2};
	}
	close OUT;
	# cdhit search
	cdhit("tmp/".$k1);
}
open OUT, ">$in.fa";
my @result = `cat tmp/*.clstr`;
foreach(@result){
	chomp;
	if(/>(.*)\.\.\.\s+\*/){
		say OUT '>',$1,"\n",$seq->{$1};
	}
}
close OUT;

system("rm -rf tmp/");

# --sub--
sub readSeq {
	my $f = shift;
	my @list = glob $f;
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
	say "$bin -i $in -d 0 -o $out -c 1.0 -n 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0";
	system("$bin -i $in -d 0 -o $out -c 1.0 -n 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0");	
}
