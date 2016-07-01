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
my %seqhash;
foreach(0..$#peparr){
	my @line = split /\t/, $peparr[$_];	
	(my $p = $line[1]) =~ s/\R|\.|\*|-//g;
	$pephash{$p} = $line[0];
	$seqhash{"seq_".$_} = $p;
}

say "# building index for chr.. seq";
my %lenSeqHash;
my $hid = '';
open FAIN, $args[0];
while(<FAIN>){
	chomp;
	if(/>(.*)/){
		$hid = $1;
	}else{
		$lenSeqHash{length($_)}{$hid} = '';
	}
}
close FAIN;

# store len -> arrindex
my $index = 0;
my @arrSeq = ();
my %arrHash = ();
foreach(sort{$a<=>$b} keys %lenSeqHash){
	$arrHash{$_} = scalar(@arrSeq) + 1 - 1;
	push @arrSeq, keys %{$lenSeqHash{$_}};
}

say '# read bed file and generate position Hash';
my %exonPosHash = ();
open IN, $bed;
while(<IN>){
	chomp;
	my @ele = split;
	(my $str = $_) =~ s/\t/_/g;
	push @{$exonPosHash{$ele[0].'_'.$ele[5].'exonLeft'}{$ele[1]+1}}, $str;
	push @{$exonPosHash{$ele[0].'_'.$ele[5].'exonRight'}{$ele[2]}}, $str;
}
close IN;

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
	# find startPos
	my ($ja,$jb,$trimaa,$orf) = (split /_/,$id)[5,6,7,8];
	$trimaa =~ s/startTrim://g;
	$orf =~ s/ORF://g;
	my ($a1, $a2, $a3) = split /,/,$ja;
	my ($b1, $b2, $b3) = split /,/,$jb;
	my $startPos = $strand eq '+' ? ($a1 + $trimaa*3 + $orf) : ($b3 - $trimaa*3 - $orf);
	if($strand eq '+' && $trimaa*3 >= ($a3 - ($a1 + $orf) + 1)){
		$startPos = $b1 + $trimaa*3 - ($a3 - ($a1 + $orf) + 1);
	}
	if($strand eq '-' && $trimaa*3 >= ($b3 - $orf - $b1 + 1)){
		$startPos = $a3 - ($trimaa*3 - ($b3 - $orf - $b1 + 1));
	}
	if($strand eq '+'){
		if($loc eq 'left'){
			if($a3 >= $startPos){
				my $s = aaindex($startPos, $a3) - 1;
				$exonloc{$id}{'0-'.$s} = $strand."\t".$loc;
			}else{
				$exonloc{$id}{'NA-NA'} = $strand."\t".$loc;	
			}
		}else{
			if($a3 >= $startPos or $startPos == $b1){
				my $s = aaindex($startPos, $a3+1) - 1;
				$exonloc{$id}{$s.'-'.(length($seq)-1)} = $strand."\t".$loc;
			}else{
				$exonloc{$id}{'NA-NA'} = $strand."\t".$loc;
				die "[ERROR] ", $_,"\t",$startPos," at line ",__LINE__,"\n";
			}
		}
	}else{
		if($loc eq 'left'){
			if($b1 <= $startPos or $startPos == $a3){
				my $s = aaindex($startPos, $b1-1) - 1;
				$exonloc{$id}{$s.'-'.(length($seq)-1)} = $strand."\t".$loc;
			}else{
				$exonloc{$id}{'NA-NA'} = $strand."\t".$loc;
                                die "[ERROR] ", $_,"\t",$startPos," at line ",__LINE__,"\n";
			}
                }else{
			if($b1 <= $startPos){
				my $s = aaindex($startPos, $b1) - 1;
				$exonloc{$id}{'0-'.$s} = $strand."\t".$loc;
			}else{
				$exonloc{$id}{'NA-NA'} = $strand."\t".$loc;
			}
                }
	}
}

say "# get peptide under fdr";
say "# peptide\tstart\tend\tAlu_HSE_exon\ttag\tid\tseq";
my %chrPepHash = ();
my %info = ();
foreach(keys %pephash){
	my $len = length($_);
	my $index = exists $arrHash{$len} ? $arrHash{$len} : 0;
	foreach my $i(@arrSeq[$index..$#arrSeq]){	
		my $index = index(uc($seq{$i}),uc($_));
		if($index != -1){
			my ($start, $end) = ($index, (length($_) + $index - 1));
			die "[ERROR] '$i' not exists in hash at line ",__LINE__,"\n" if !exists $exonloc{$i};
			my $k = join ";", keys %{$exonloc{$i}};
			my $tag = 'N';
			foreach my $k2(keys %{$exonloc{$i}}){
				next if $k2 =~ /NA/;
				my ($s, $e) = split /-/,$k2;
				my $min = min($s, $e, $start, $end);
				my $max = max($s, $e, $start, $end);
				$tag = 'Y' if abs($end - $start) + 1 + abs($e - $s) + 1 > ($max - $min + 1);
			}
			$chrPepHash{$_} = $pephash{$_};
			my ($chr, $jl, $jr, $strand, $annol, $annor, $startTrim, $ORF) = (split /_/,$i)[0,1,2,4,5,6,7,8];
			$startTrim =~ s/startTrim://g;
			$ORF =~ s/ORF://g;
			$startTrim += $start;
			my $l = (split /,/, $annol)[0];
			my $r = (split /,/, $annor)[2];
			my $startPos = $strand eq '+' ? ($l+$startTrim*3+$ORF) : ($r-$startTrim*3-$ORF);
			if($strand eq '+' && $startPos > $jl){
                        	$startPos = $jr + $startPos - $jl - 1;              
			}
			if($strand eq '-' && $startPos < $jr){
                                $startPos = $jl - ($jr - $startPos) + 1;
                        }
			my $endPos = $strand eq '+' ? ($startPos-1+length($_)*3) : ($startPos+1-length($_)*3);
			if($strand eq '+' && $startPos <= $jl && $endPos > $jl){
				$endPos = $jr + $endPos - $jl - 1;
			}
			if($strand eq '-' && $startPos >= $jr && $endPos < $jr){
                                $endPos = $jl - ($jr - $endPos) + 1;
                        }
			#say $_,"\t",$startTrim,"\t",$ORF,"\t",$startPos,"\t",$endPos,"\t",$i if $_ =~ /RGLSQSALPYRR/;
			my $exons = '';
			if(exists $exonPosHash{$chr.'_'.$strand.'exonLeft'}{$jr}){
				my @arr = @{$exonPosHash{$chr.'_'.$strand.'exonLeft'}{$jr}};
				say 'right: ',join ';',@arr if $_ =~ /RGLSQSALPYRR/;
				foreach my $a(@arr){
					my ($minLeft,$maxRight) = (split /_/,$a)[1,2];
					$minLeft = $minLeft + 1;
					my @p = sort($maxRight, $minLeft, $startPos, $endPos);
					my $dis = abs($maxRight-$minLeft)+1+ abs($endPos-$startPos)+1;
					my $overlapStatus = abs($p[3]-$p[0])+1 < $dis ? 1 : 0;
					$exons = $a.';'.$exons if $maxRight >= $startPos && $maxRight >= $endPos && $overlapStatus == 1;
				}
			}
			if(exists $exonPosHash{$chr.'_'.$strand.'exonRight'}{$jl}){
				my @arr = @{$exonPosHash{$chr.'_'.$strand.'exonRight'}{$jl}};
				say 'left: ',join ';',@arr if $_ =~ /RGLSQSALPYRR/;
                                foreach my $a(@arr){
                                        my ($minLeft,$maxRight) = (split /_/,$a)[1,2];
					$minLeft = $minLeft + 1;
					my @p = sort($maxRight, $minLeft, $startPos, $endPos);
                                        my $dis = abs($maxRight-$minLeft)+1+ abs($endPos-$startPos)+1;
                                        my $overlapStatus = abs($p[3]-$p[0])+1 < $dis ? 1 : 0;
					$exons = $a.';'.$exons if $minLeft <= $startPos && $minLeft <= $endPos && $overlapStatus == 1;
                                }
                        }
			push @{$info{$_}}, $_."\t".$start."\t".$end."\t".$k."\t".$tag."\t".$i."\t".$seq{$i}."\t".$exons."\t".$startPos."\t".$endPos;
		}
	}
}

foreach my $peptide(keys %info){
	my @k = @{$info{$peptide}};
	#say $_,"\t",$chrPepHash{$peptide} for @k;
}

my $localResult = localFDR(\%chrPepHash);
my @arr = split /\R/, $localResult;
foreach(@arr){
	my @ele = split /\t/;
	my @k = @{$info{$ele[0]}};
	say $_,"\t",$chrPepHash{$ele[0]} for @k;
}


#system("rm -rf tmp/");

sub aaindex {
	my ($start, $end) = @_;
	my $len = abs($end - $start) + 1;
	my $index = -1;
	if($len % 3 == 0 ){
		$index = $len / 3;
	}else{
		$index = ($len - ($len%3))/3 + 1;
	}
	return $index;	
}

sub usage {
my $usage = "USAGE: ".basename($0)." junctionpep_path percolator_path bedfile_path
Example:
".basename($0)." result_v5/Alu/junctionPep/GM19207.rna.fa ~/comet/GM19207_tsv/0.05/ data/Ensembl_Alu_25bp_0.5_CDS.unique.sorted.bed\n";
return $usage;
}

sub localFDR {
	my $chrPepHash = shift;
	my @uniq = sort {$a<=>$b} uniq(values %{$chrPepHash});
	foreach my $k(keys %{$chrPepHash}){
		#say $k,"\t",$chrPepHash->{$k};
	}	
	my $result = '';
	my $FDR = 0.01;
	foreach my $s (@uniq){
		my $sum = 0;
		my $content = '';
		my $num = 0;
		my $total = scalar(keys %{$chrPepHash});
		foreach my $l (keys %{$chrPepHash}){
			my $pep = $chrPepHash->{$l};
			if($pep <= $s){
	                        $sum += $pep;
				$num += 1;
                        	$content .= $l."\t".$pep."\n";
                	}
		}
		my $fdr = $sum / $num;
		last if $fdr > $FDR;
		$result = $content;
	}
	return $result;
}
