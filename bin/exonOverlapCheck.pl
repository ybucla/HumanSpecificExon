#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use List::Util qw( min max);
use List::MoreUtils qw(uniq);
use 5.010;

my $resultfile = '../GM18486.0.05.list';# '../result_v5/LCLs/Alu_all.0.05_0.01.list';
my $exonbedfile = '../data/Ensembl/exon.unique.bed';

# result 2 bed
my $hash = result2bed($resultfile);

# get overlaped exon with same junction
foreach(keys %{$hash}){
	next if !/GM18486/;
	my $junctionbed = 'tmp/tmp_'.$_;
	my @r = `bedtools intersect -a $junctionbed -b $exonbedfile -s -wo`;
	my %result = ();
	foreach my $line(@r){
		my @ele = split /\t/, $line;
		my $junctionLen = $ele[2] - $ele[1];
		my $overlapLen = $ele[12];
		my $id = join "_",(split /_/,$ele[3])[0..7];
		my $pep = (split /_/,$ele[3])[9];
		next if $overlapLen != $junctionLen;
		$result{$_}{$pep}{$id}{(join ":",@ele[6,9,11])} = '';
	}
	# check whether a pep could be mapped to mutiple exon peptides? if so discard it, else keep
	foreach my $pep(keys %{$result{$_}}){
		foreach my $id(keys %{$result{$_}{$pep}}){
			my $len = scalar(keys %{$result{$_}{$pep}{$id}});
			#next if $len > 1;
			say $hash->{$_}{$pep}{$id},"\t",$len;
		}
	}
}


sub result2bed {
	my $in = shift; # 

	my $head = '';
	my %hash = ();
	open IN, $in;
	while(<IN>){
		chomp;
		($head = $_) =~ s/======|\.rna//g if /==/;
		next if /#|=/;
		my ($pep, $s, $e, $splice, $tag, $id, $seq) = split;
		next if $tag eq 'N';
		$hash{$head}{$pep}{$id} = $_;
	}
	close IN;
	# write junction (Y) to bed
	my $tmpdir = 'tmp';
	mkdir $tmpdir if !-e $tmpdir;
	foreach(keys %hash){
		my $tmpout = $tmpdir.'/tmp_'.$_;
		open OUT, '>'.$tmpout;
		foreach my $pep(keys %{$hash{$_}}){
			foreach my $i(keys %{$hash{$_}{$pep}}){
				my $line = $hash{$_}{$pep}{$i};
				my ($peptide, $s, $e, $splice, $tag, $id, $seq) = split /\t/, $line;	
				my @arr = split /_/,$id;
				my ($chr, $strand) = @arr[0,4];
				my ($left, $jl) = (split /,/,$arr[5])[0,2];
				my ($jr, $right) = (split /,/,$arr[6])[0,2];
				foreach my $ss((split/;/,$splice)){
					my ($splice_s, $splice_e) = split /-/, $ss;
					my @num = sort($s, $e, $splice_s, $splice_e);
					my $p = substr($seq, $num[1], $num[2]-$num[1]+1);
					if($splice_s > 0){
						if($strand eq '+'){
							say OUT join "\t", ($chr,$jr-1,$right,$id."_right_".$peptide."_".$_,0,$strand);
						}else{
							say OUT join "\t", ($chr,$left-1,$jl,$id."_right_".$peptide."_".$_,0,$strand);
						}
					}else{
						if($strand eq '+'){
							say OUT join "\t", ($chr,$left-1,$jl,$id."_left_".$peptide."_".$_,0,$strand);
						}else{
							say OUT join "\t", ($chr,$jr-1,$right,$id."_left_".$peptide."_".$_,0,$strand);
						}
					}
				}
			}
		}
		close OUT;
	}
	return \%hash;
}

sub translate {
	my $seq = shift;
	$seq = uc($seq);
	my %aacode = (
		TTT => "F", TTC => "F", TTA => "L", TTG => "L",
		TCT => "S", TCC => "S", TCA => "S", TCG => "S",
		TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
		TGT => "C", TGC => "C", TGA => "*", TGG => "W",
		CTT => "L", CTC => "L", CTA => "L", CTG => "L",
		CCT => "P", CCC => "P", CCA => "P", CCG => "P",
		CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
		CGT => "R", CGC => "R", CGA => "R", CGG => "R",
		ATT => "I", ATC => "I", ATA => "I", ATG => "M",
		ACT => "T", ACC => "T", ACA => "T", ACG => "T",
		AAT => "N", AAC => "N", AAA => "K", AAG => "K",
		AGT => "S", AGC => "S", AGA => "R", AGG => "R",
		GTT => "V", GTC => "V", GTA => "V", GTG => "V",
		GCT => "A", GCC => "A", GCA => "A", GCG => "A",
		GAT => "D", GAC => "D", GAA => "E", GAG => "E",
		GGT => "G", GGC => "G", GGA => "G", GGG => "G",
	); # this is the hash table for the amino acids
	my @result = ();
	my $len = length($seq);
	for(my $orf = 1;$orf<4;$orf++){
		my $pep = '';
		for(my $i=$orf-1; $i<$len-2;$i=$i+3){
			$pep .= $aacode{substr($seq, $i, 3)};
		}
		push @result, $pep;
	}
	return \@result;
}
