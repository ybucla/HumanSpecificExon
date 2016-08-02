#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use List::Util qw( min max);
use List::MoreUtils qw(uniq);
use 5.010;
use Thread;

my $start_run = time();

# -----------------------------
my @args = @ARGV or die usage(); # result_v5/Alu/junctionPep/GM18486.rna.fa result_v5/Alu/cometout/GM18486.rna/0.05
die usage() if scalar(@args) < 2;
my $bed = $args[2]; # "/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/Ensembl_Alu_25bp_0.5_CDS.unique.sorted.bed";
my $SJdir = $args[3]; # /u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/SJ_out/TISSUE
my $threadnum = 1;

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

my %junctionHash = ();
say '# building index for junction SJ.out.tab';
my @SJfiles = glob "$SJdir/$sample*";
foreach my $f(@SJfiles){
	open IN, $f;
	while(<IN>){
		chomp;
		my @ele = split /\t/;
		my $strand = $ele[3] eq '1' ? '+' : '-';
		my $junction = $ele[0].'_'.($ele[1]-1).'_'.($ele[2]+1).'_'.$strand;
		$junctionHash{$junction} = $ele[5];
	}
	close IN;
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
	my $seq = $seq{$id};
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
				my $s = $startPos == $b1 ? 0 : (aaindex($startPos, $a3+1) - 1);
				$exonloc{$id}{$s.'-'.(length($seq)-1)} = $strand."\t".$loc;
			}else{
				$exonloc{$id}{'NA-NA'} = $strand."\t".$loc;
				die "[ERROR] ", $_,"\t",$startPos," at line ",__LINE__,"\n";
			}
		}
	}else{
		if($loc eq 'left'){
			if($b1 <= $startPos or $startPos == $a3){
				my $s = $startPos == $a3 ? 0 : (aaindex($startPos, $b1-1) - 1);
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


say "# build index2 for chr.. seq";
my %tmplen = ();
my %tmpfirst = ();
$tmpfirst{substr($_,0,1)}='' for keys %pephash;
foreach(keys %seq){
	my $s = uc($seq{$_});
	foreach my $f (keys %tmpfirst){
		my $i = index($s, $f);
		next if $i == -1;
		my $l = length($s)- $i;
		$tmplen{$f}{$l}{$_} = '';
	}
}
my %arrSeq2 = ();
my %arrStart = ();
foreach my $f(keys %tmplen){
	foreach my $l(sort{$a<=>$b} keys %{$tmplen{$f}}){
		if(exists $arrSeq2{$f}){
			$arrStart{$f}{$l} = scalar(@{$arrSeq2{$f}}) + 1 - 1;
			push @{$arrSeq2{$f}}, keys %{$tmplen{$f}{$l}};
		}else{
			$arrStart{$f}{$l} = 0;
			my @t = keys %{$tmplen{$f}{$l}};
			$arrSeq2{$f} = \@t;
		}
	}
}

# arr index
my $indexThread :shared = 0;
# result hash
my %chrPepHash :shared = ();
my @infoArr :shared = ();
my %info = ();

say "# multiple thread to get junction peptides";
my @threads;
my $tempcount = 0;
foreach(1..$threadnum){
	$threads[$tempcount]=Thread->new(\&run,\%pephash, \%arrHash, \%arrStart, \%arrSeq2, \%seq, \%exonloc, \%exonPosHash, $_);
	$tempcount++;
}

my $r = $_->join() for @threads;

foreach(@infoArr){
	my @ele = split /\t/,$_;
	push @{$info{$ele[0]}}, $_;
}

say "# get peptide under fdr";
say "# peptide\tstart\tend\tAlu_HSE_exon\ttag\tid\tseq";
my %novelChrPepHash = ();
foreach my $peptide(keys %info){
	my @k = @{$info{$peptide}};
	foreach(@k){
		my $head = (split /\t/,$_)[5];
		my @arr = split /_/, $head;
		my $junc = join '_', @arr[0,1,2,4];
		die "ERROR:\t junction not exists in SJ.out.tab file\n$junc\n" if !exists $junctionHash{$junc};
		if($junctionHash{$junc} == 0){
			$novelChrPepHash{$peptide} = $chrPepHash{$peptide};
		}
	}
	#say $_,"\t",$chrPepHash{$peptide} for @k;
}
warn "# num of pep from novel junction:\t",scalar(keys %novelChrPepHash),"\n";
warn "# num of pep from annot junction:\t", scalar(keys %chrPepHash)-scalar(keys %novelChrPepHash),"\n";

foreach my $peptide(keys %chrPepHash){
	next if exists $novelChrPepHash{$peptide};
	my @k = @{$info{$peptide}};
	say $_,"\t",$chrPepHash{$peptide},"\t1" for @k;
}

my $localResult = localFDR(\%novelChrPepHash);
my @arr = split /\R/, $localResult;
foreach(@arr){
	my @ele = split /\t/;
	my @k = @{$info{$ele[0]}};
	say $_,"\t",$chrPepHash{$ele[0]},"\t0" for @k;
}


#system("rm -rf tmp/");

# ---------------------------------------
my $end_run = time();
my $run_time = $end_run - $start_run;
warn "Job took $run_time seconds\n";

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
	warn "# glob FDR 5%:\t",scalar(keys %{$chrPepHash}),"\n";
	foreach my $k(keys %{$chrPepHash}){
		#warn $k,"\t",$chrPepHash->{$k},"\n";
	}	
	my $result = '';
	my $FDR = 0.01;
	my $n = 0;
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
		$n++;
		$result = $content;
	}
	warn "# local FDR 1%:\t", $n,"\n";
	return $result;
}

sub run {
my ($pephash, $arrHash, $arrStart, $arrSeq2, $seq, $exonloc, $exonPosHash, $threadNum )= @_;

my @arrPepHash = sort {$a cmp $b} keys %{$pephash};
while( (my $j = getindex()) < scalar(@arrPepHash)){
	my $pepj = $arrPepHash[$j];
	my $len = length($pepj);
	my $index = exists $arrHash->{$len} ? $arrHash->{$len} : 0;
	my $firstAA = substr($pepj,0,1);
	my $iStart = exists $arrStart->{$firstAA}{$len} ? $arrStart->{$firstAA}{$len} : 0;
	my $iEnd = scalar(@{$arrSeq2->{$firstAA}})-1;
	warn "Thread-$threadNum: ",$j,"\t",$len,"\t",$firstAA,"\t",$iStart,"\t",$iEnd,"\n";
	my @arr = @{$arrSeq2->{$firstAA}};
	foreach my $i(@arr[$iStart..$iEnd]){
	#foreach my $i(@arrSeq[$index..$#arrSeq]){
		die "ERROR:\t",$j,"\t",$pepj,"\n" if !defined($pepj);
		my $index = index(uc($seq->{$i}),uc($pepj));
		if($index != -1){
			my ($start, $end) = ($index, (length($pepj) + $index - 1));
			die "[ERROR] '$i' not exists in hash at line ",__LINE__,"\n" if !exists $exonloc->{$i};
			my $k = join ";", keys %{$exonloc->{$i}};
			my $tag = 'N';
			foreach my $k2(keys %{$exonloc->{$i}}){
				next if $k2 =~ /NA/;
				my ($s, $e) = split /-/,$k2;
				my $min = min($s, $e, $start, $end);
				my $max = max($s, $e, $start, $end);
				$tag = 'Y' if abs($end - $start) + 1 + abs($e - $s) + 1 > ($max - $min + 1);
			}
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
			my $endPos = $strand eq '+' ? ($startPos-1+length($pepj)*3) : ($startPos+1-length($pepj)*3);
			if($strand eq '+' && $startPos <= $jl && $endPos > $jl){
				$endPos = $jr + $endPos - $jl - 1;
			}
			if($strand eq '-' && $startPos >= $jr && $endPos < $jr){
                                $endPos = $jl - ($jr - $endPos) + 1;
                        }
			my $exons = '';
			if(exists $exonPosHash->{$chr.'_'.$strand.'exonLeft'}{$jr}){
				my @arr = @{$exonPosHash->{$chr.'_'.$strand.'exonLeft'}{$jr}};
				foreach my $a(@arr){
					my ($minLeft,$maxRight) = (split /_/,$a)[1,2];
					$minLeft = $minLeft + 1;
					my @p = sort($maxRight, $minLeft, $startPos, $endPos);
					my $dis = abs($maxRight-$minLeft)+1+ abs($endPos-$startPos)+1;
					my $overlapStatus = abs($p[3]-$p[0])+1 < $dis ? 1 : 0;
					$exons = $a.';'.$exons if $maxRight >= $startPos && $maxRight >= $endPos && $overlapStatus == 1;
				}
			}
			if(exists $exonPosHash->{$chr.'_'.$strand.'exonRight'}{$jl}){
				my @arr = @{$exonPosHash->{$chr.'_'.$strand.'exonRight'}{$jl}};
                                foreach my $a(@arr){
                                        my ($minLeft,$maxRight) = (split /_/,$a)[1,2];
					$minLeft = $minLeft + 1;
					my @p = sort($maxRight, $minLeft, $startPos, $endPos);
                                        my $dis = abs($maxRight-$minLeft)+1+ abs($endPos-$startPos)+1;
                                        my $overlapStatus = abs($p[3]-$p[0])+1 < $dis ? 1 : 0;
					$exons = $a.';'.$exons if $minLeft <= $startPos && $minLeft <= $endPos && $overlapStatus == 1;
                                }
                        }
			{
				lock(%chrPepHash);
				$chrPepHash{$pepj} = $pephash->{$pepj};
			}
			{
				lock(@infoArr);
				my $str = $pepj."\t".$start."\t".$end."\t".$k."\t".$tag."\t".$i."\t".$seq->{$i}."\t".$exons."\t".$startPos."\t".$endPos;
				push @infoArr, $str;
			}
		}
	}
}
}

sub getindex {

	lock($indexThread);

	return $indexThread++;
}
