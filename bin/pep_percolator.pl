#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use List::Util qw( min max);
use List::MoreUtils qw(uniq);
use 5.010;

my @args = @ARGV or die usage(); # result_v5/Alu/junctionPep/GM18486.rna.fa result_v5/Alu/cometout/GM18486.rna/0.05
die usage() if scalar(@args) < 2;

my $bed = $args[2]; # "/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/Ensembl_Alu_25bp_0.5_CDS.unique.sorted.bed";

(my $sample = basename($args[0])) =~ s/\.fa//g;;
my $fdr = basename($args[1]);

say "======$sample";
say "# obtain percolator result from '$sample' with fdr '$fdr'";
say "# bedfile: $bed";
my @peparr = `awk -F '\t' '\$3 < $fdr' $args[1]/*.pep.percolator | cut -f 4,5`;
my %pephash;
foreach(0..$#peparr){
	my @line = split /\t/, $peparr[$_];	
	(my $p = $line[1]) =~ s/\R|\.|\*|-//g;
	$pephash{$p} = $line[0];
}

say "# read junction pep and write to bed";
my %seq;
my $head = '';
mkdir "tmp" if !-e "tmp";
open BED,">tmp/$sample.bed";
open IN, $args[0];
while(<IN>){
	chomp;
	if(/>(.*)/){
		$head = $1;
		my @arr = split /_/,$head;
		my ($chr, $strand) = @arr[0,4];	
		my ($left, $jl) = (split /,/,$arr[5])[1,2];
		my ($right, $jr) = (split /,/,$arr[6])[1,0];
		say BED $chr,"\t",$left-1,"\t",$jl,"\t",$head,"_left","\t0\t",$strand;
		say BED $chr,"\t",$jr-1,"\t",$right,"\t",$head,"_right","\t0\t",$strand;
	}else{
		$seq{$head} = $_;
	}
}
close IN;
close BED;

say "# bedtools to find Exon location";
my %exonloc;
my @r = `bedtools coverage -a tmp/$sample.bed -b $bed -s`;
foreach(@r){
	chomp;
	my @line = split;
	next if $line[7] == 0;
	my $strand = $line[5];
	my $loc = (split /_/, $line[3])[-1];
	(my $id = $line[3]) =~ s/_left|_right//g;
	(my $seq = $seq{$id}) =~ m/[a-z]/g;
	my $splice = pos($seq) - 1;
	if($strand eq '+'){
		if($loc eq 'left'){			
			$exonloc{$id}{'0-'.$splice} = $strand."\t".$loc;
		}else{
			$exonloc{$id}{$splice.'-'.(length($seq)-1)} = $strand."\t".$loc;
		}
	}else{
		if($loc eq 'left'){
			$exonloc{$id}{$splice.'-'.(length($seq)-1)} = $strand."\t".$loc;
                }else{
			$exonloc{$id}{'0-'.$splice} = $strand."\t".$loc;
                }
	}
}

#foreach my $k1(keys %exonloc){
#	foreach my $k2(keys %{$exonloc{$k1}}){
#		say $k1, "\t", $k2, "\t", $exonloc{$k1}{$k2};
#	}
#}

say "# get peptide under fdr";
say "# peptide\tstart\tend\tsplice\ttag\tid\tseq";
my %chrPepHash = ();
my %info = ();
foreach(keys %pephash){
	foreach my $i(keys %seq){
		my $index = index(uc($seq{$i}),uc($_));
		if($index != -1){
			my ($start, $end) = ($index, (length($_) + $index - 1));
			my $k = join ";", keys %{$exonloc{$i}};
			my $tag = 'N';
			foreach my $k2(keys %{$exonloc{$i}}){
				my ($s, $e) = split /-/,$k2;
				my $min = min($s, $e, $start, $end);
				my $max = max($s, $e, $start, $end);
				$tag = 'Y' if ($end - $start + 1 + $e - $s + 1) > ($max - $min + 2);
			}
			#say $_,"\t",$start,"\t",$end,"\t",$k,"\t",$tag,"\t",$i,"\t",$seq{$i};
			$chrPepHash{$_} = $pephash{$_};
			push @{$info{$_}}, $_."\t".$start."\t".$end."\t".$k."\t".$tag."\t".$i."\t".$seq{$i};
		}
	}
}

foreach my $peptide(keys %info){
	my @k = @{$info{$peptide}};
	say $_,"\t",$chrPepHash{$peptide} for @k;
}

my $localResult = localFDR(\%chrPepHash);
my @arr = split /\R/, $localResult;
foreach(@arr){
	my @ele = split /\t/;
	my @k = @{$info{$ele[0]}};
	#say $_,"\t",$chrPepHash{$ele[0]} for @k;
}


#system("rm -rf tmp/");

sub usage {
my $usage = "USAGE: ".basename($0)." junctionpep_path percolator_path bedfile_path
Example:
".basename($0)." result_v5/Alu/junctionPep/GM19207.rna.fa ~/comet/GM19207_tsv/0.05/ data/Ensembl_Alu_25bp_0.5_CDS.unique.sorted.bed\n";
return $usage;
}

sub localFDR {
	my $chrPepHash = shift;
	my @uniq = sort {$a<=>$b} uniq(values %{$chrPepHash});
	my $result = '';
	my $FDR = 0.01;
	foreach my $s (@uniq){
		my $sum = 0;
		my $content = '';
		my $total = scalar(keys %{$chrPepHash});
		foreach my $l (keys %{$chrPepHash}){
			my $pep = $chrPepHash->{$l};
			if($pep < $s){
	                        $sum += $pep;
                        	$content .= $l."\t".$pep."\n";
                	}
		}
		my $fdr = $sum / $total;
		last if $fdr > $FDR;
		$result = $content;
	}
	return $result;
}
