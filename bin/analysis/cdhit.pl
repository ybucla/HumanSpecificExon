#! /u/home/y/ybwang/perl
use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my $in = shift or die "ERROR: $0 lymph.0.05.list fasta\nPlease input '0.05.list file'\n"; #'result_v5/TISSUE/HSE/nolocalFDR/lymph.nolocalFDR.0.05.list';
my $fasta = shift or die "ERROR: $0 lymph.0.05.list fasta\nPlease input 'fasta' file\n" ;#"result_v5/TISSUE/HSE/junctionPep/lymph.fa";

die "ERROR: '$in' file not exists\n" if !-e $in;
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
	next if /^#|^===|^\//;
	my ($tag, $j) = (split /\t/)[4,5];
	next if $tag eq 'N';	
	my @arr = split /,|_/,$j;
	my $id = join "_", @arr[0,1,2,4];	
	$junction{$id}{$j."_".$head} = $seq->{$j."_".$head};
	#$junction{$id}{$j} = $seq->{$j};
}
close IN;


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

# --sub--
sub readSeq {
	my @list = glob $fasta;
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
				#$seq{$id} = $_;
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
	say "$bin -i $in -d 0 -o $out -c 0.8 -n 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0";
	system("$bin -i $in -d 0 -o $out -c 0.8 -n 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0");	
}
