#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use List::Util qw( min max);
use List::MoreUtils qw(uniq);
use 5.010;

my $resultfile = shift or die USAGE();# 'result_v5/LCLs/Alu_all/Alu_all.0.05.list';# '../result_v5/LCLs/Alu_all.0.05_0.01.list';
my $exonbedfile = './data/Ensembl/exon.unique.bed';

die USAGE() if !-e $resultfile;

# result 2 bed
my $hash = result2bed($resultfile);

# get overlaped exon with same junction
foreach(sort{$a cmp $b} keys %{$hash}){
	#next if !/GM18870/;
	say '======'.$_;
	say "# obtain exon overlap result from '$_'";
	say "# bedfile: $exonbedfile";
	say "# peptide\tstart\tend\tsplice\ttag\tid\tseq\tqvalue\toverlap";
	my $junctionbed = 'tmp/tmp_'.$_;
	my @r = `bedtools intersect -a $junctionbed -b $exonbedfile -s -wo`;
	
	my %result = ();
	foreach my $line(@r){
		chomp($line);
		my @ele = split /\t/, $line;
		my $junctionLen = $ele[2] - $ele[1];
		my $overlapLen = $ele[12];
		my $id = join "_",(split /_/,$ele[3])[0..8];
		my $strand = $ele[5];
		my $pep = (split /_/,$ele[3])[11];
		next if $overlapLen != $junctionLen;
		if($pep =~ /[a-z]/){ # if the comet peptide cross the junction sites
			if(($ele[3] =~ /left/ && $strand eq '+') or ($ele[3] =~ /right/ && $strand eq '-')){
				$result{$_}{$pep}{$id}{(join ":",@ele[6,9,11])} = '' if $ele[2] == $ele[8];
			}
			if(($ele[3] =~ /left/ && $strand eq '-') or ($ele[3] =~ /right/ && $strand eq '+')){
				$result{$_}{$pep}{$id}{(join ":",@ele[6,9,11])} = '' if $ele[1] == $ele[7];
			}
		}else{ # if the comet peptide fall into one side of junction
			$result{$_}{$pep}{$id}{(join ":",@ele[6,9,11])} = '';
		}
	}
	# check whether a pep could be mapped to mutiple exon peptides? if so discard it, else keep
	foreach my $pep(keys %{$result{$_}}){
		foreach my $id(keys %{$result{$_}{$pep}}){
			my $len = scalar(keys %{$result{$_}{$pep}{$id}});
			next if $len > 1;
			say $hash->{$_}{$pep}{$id},"\t",$len;
		}
	}
}

#system("rm -rf tmp/");

sub result2bed {
	my $in = shift; # 

	my $head = '';
	my %hash = ();
	open IN, $in;
	while(<IN>){
		chomp;
		($head = $_) =~ s/======|\..*$//g if /==/;
		next if /#|=/;
		my ($pep, $s, $e, $splice, $tag, $id, $seq, $qvalue) = split;
		next if $tag eq 'N';
		$pep = substr($seq, $s, $e-$s+1);
		$hash{$head}{$pep}{$id} = join "\t", ($pep, $s, $e, $splice, $tag, $id, $seq, $qvalue);
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
				my ($startTrim, $ORF) = @arr[7,8];
				$startTrim =~ s/startTrim://g;
				$ORF =~ s/ORF://g;
				$startTrim = $startTrim + $s;
				my $endPos = 1;
				foreach my $ss((split/;/,$splice)){
					my ($splice_s, $splice_e) = split /-/, $ss;
					next if $e < $splice_s or $splice_e < $s;
					my @num = sort($s, $e, $splice_s, $splice_e);
					if($pep =~ /[a-z]/){ # peptide include splicing site
						my $startPos = $strand eq '+' ? ($left+$ORF+$startTrim*3) : ($right-$ORF-$startTrim*3);
						if($strand eq '+'){
							$endPos = length($peptide) * 3 - ($jl - $startPos + 1) + $jr - 1;
							if($splice_s > 0){
								say OUT join "\t", ($chr,$jr-1,$endPos,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
							}else{
								say OUT join "\t", ($chr,$startPos-1,$jl,$id."_index:$s"."_left_".$pep."_".$_,0,$strand);
							}
						}else{							
							$endPos = $jl + 1 - (length($peptide) * 3 - ($startPos - $jr + 1));
							if($splice_s > 0){
								say OUT join "\t", ($chr,$endPos-1,$jl,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
							}else{
								say OUT join "\t", ($chr,$jr-1,$startPos,$id."_index:$s"."_left_".$pep."_".$_,0,$strand);
							}
						}
					}else{ # peptide include no splicing site
						if($strand eq '+'){
							my $startPos = $left+$ORF+$startTrim*3;
							if($startTrim*3 - 1 >= abs($jl - $left - $ORF)){
								$startPos = $jr + $startTrim*3 - abs( $jl - $left - $ORF) - 1 - 1 + 1;
							}			
							$endPos = $startPos - 1 + length($peptide) * 3;
							if($splice_s > 0){
								say OUT join "\t", ($chr,$startPos-1,$endPos,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
							}else{
								say OUT join "\t", ($chr,$startPos-1,$endPos,$id."_index:$s"."_left_".$pep."_".$_,0,$strand);
							}
						}else{
							my $startPos = $right-$ORF-$startTrim*3;
							if($startTrim*3 - 1 >= abs($right - $ORF - $jr)){
								$startPos = $jl - ($startTrim*3 - abs($right - $jr - $ORF) - 1) + 1 - 1;
							}
							$endPos = $startPos + 1 - length($peptide) * 3;
							if($splice_s > 0){
								say OUT join "\t", ($chr,$endPos-1,$startPos,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
							}else{
								say OUT join "\t", ($chr,$endPos-1,$startPos,$id."_index:$s"."_left_".$pep."_".$_,0,$strand);
							}
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

sub USAGE {
my $str = "USAGE: $0 file.list(bin/tjScript result)\n";
return $str;
}
