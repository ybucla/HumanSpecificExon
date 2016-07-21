#! /u/home/y/ybwang/perl

# 

use strict;
use warnings;
use File::Basename;
use Data::Table;
use List::Util qw( min max );
use List::MoreUtils qw(uniq);
use 5.010;

my $resultfile = shift or die USAGE();
my $bedfile = shift or die USAGE();
my $SJfile = shift || 'data/SJ_out/GM18486.rna.SJ';
my $Aluelement = shift || 'false';

die "[ERROR] '$resultfile' not exists!\n".USAGE() if !-e $resultfile;
die "[ERROR] '$bedfile' not exists!\n".USAGE() if !-e $bedfile;

# file build-in
my $exonbedfile = './data/Ensembl/exon.unique.bed';
# remove 'tmp' dir?
my $rmtmp = 0;

# print parameters
warn "[PARAMETERS] list_file:\t", $resultfile,"\n";
warn "[PARAMETERS] bed_file:\t", $bedfile,"\n";
warn "[PARAMETERS] SJ_file:\t", $SJfile,"\n";
warn "[PARAMETERS] Alu_element:\t", $Aluelement,"\n";

# read SJ file
warn "[log]# reading SJ.out.file\n";
my $SJhash = SJ2hash($SJfile);
# read bed file
warn "[log]# reading bedfile\n";
my $bedhash = bed2hash($bedfile);
# result 2 bed
warn "[log]# reading list file\n";
my $hash = result2bed($resultfile);
# cat tow bed file
mkdir 'tmp' if !-e 'tmp';
my @cat = `cat $exonbedfile $bedfile > tmp/exon.bed`;


# get overlaped exon with same junction
warn "[log]# analysis peptide and exon overlaping\n";
foreach(sort{$a cmp $b} keys %{$hash}){	
	#next if !/GM18486/;
	warn "[log]# ======$_\n";
	say '======'.$_;
	say "# obtain exon overlap result from '$_'";
	say "# bedfile: $exonbedfile";
	say "# peptideChr\tstart\tend\tid\t\tscore\tstrand\toverlapexon\ttargetExon\toverlapexonNum";
	my $junctionbed = 'tmp/tmp_'.$_;
	# findout Alu element overlaped peptide
	my %peptideTypeHash;
	my @Alu_overlap = `bedtools coverage -a $junctionbed -b data/hg19_repeatMasker_Alu.sorted_pos.bed`;
	foreach my $line(@Alu_overlap){
		chomp($line);
		my @ele = split /\t/,$line;
		my $str = join "\t",@ele[0..5];
		$ele[9] > 0 ? ($peptideTypeHash{$str} = "Aluelement") : ($peptideTypeHash{$str} = "Aluexon");
		$peptideTypeHash{$str} = "Aluexon" if $Aluelement ne 'true';
	}
	
	my @r = `bedtools intersect -a $junctionbed -b tmp/exon.bed -s -wo`;	
	my %result = ();
	my %uni = ();
	foreach my $line(@r){
		chomp($line);
		my @ele = split /\t/, $line;
		my $junctionLen = $ele[2] - $ele[1];
		my $overlapLen = $ele[12];
		next if $overlapLen != $junctionLen;
		my $id = $ele[3];#join "_",(split /_/,$ele[3])[0..8];
		my $strand = $ele[5];
		my $pep = (split /_/,$ele[3])[11];
		my $pepCor = join "\t",(@ele[0..2],$id,@ele[4,5]);
		my $str = (join "_",@ele[6,7,8,9,10,11]);
		my $exonstr = join "\t",@ele[0..5];
		# find cross start
		my @arr = split /_/,$id;
		my ($left, $jl) = (split /,/,$arr[5])[0,2];
                my ($jr, $right) = (split /,/,$arr[6])[0,2];
                my ($startTrim, $ORF, $index) = @arr[7,8,9];
                $startTrim =~ s/startTrim://g;
                $index =~ s/index://g;
                $ORF =~ s/ORF://g;
		$startTrim = $startTrim + $index;
		my $crossSplice = "false";
                my ($a1, $a,  $a3) = ($left, 0, $jl);
                my ($b1, $b2, $b3) = ($jr, 0, $right);
                my $startPos = $strand eq '+' ? ($a1+$startTrim*3+$ORF) : ($b3-$startTrim*3-$ORF);
                if($strand eq '+' && $startTrim*3 >= ($a3 - ($a1 + $ORF) + 1)){
                        $startPos = $b1 + $startTrim*3 - ($a3 - ($a1 + $ORF) + 1);              
                }
                if($strand eq '-' && $startTrim*3 >= ($b3 - $ORF - $b1 + 1)){
                	$startPos = $a3 - ($startTrim*3 - ($b3 - $ORF - $b1 + 1));
                }
		my $endTag = $strand eq '+' ? ($startPos+length($pep)*3-1) : ($startPos-length($pep)*3+1);
                $crossSplice = 'true' if $strand eq '+' && $startPos <= $a3 && $endTag > $a3;
		$crossSplice = 'true' if $strand eq '-' && $startPos >= $b1 && $endTag < $b1;
		# find cross end
		if($crossSplice eq 'true'){ # if the comet peptide cross the junction sites
                        if(($ele[3]=~/left/ && $strand eq '+') or ($ele[3]=~/right/ && $strand eq '-')){
				if($ele[2] == $ele[8]){
					if($peptideTypeHash{$exonstr} eq "Aluelement"){
        		                	$uni{$_}{$pepCor}{$str} = '';
			                }else{
                	        		$uni{$_}{$pepCor}{$str} = '' if $ele[9] !~ /exon/;
		        	        }
				}
                        }
                        if(($ele[3]=~/left/ && $strand eq '-') or ($ele[3]=~/right/ && $strand eq '+')){
				if($ele[1] == $ele[7]){
					if($peptideTypeHash{$exonstr} eq "Aluelement"){
                                                $uni{$_}{$pepCor}{$str} = '';
                                        }else{
                                                $uni{$_}{$pepCor}{$str} = '' if $ele[9] !~ /exon/;
                                        }
				}
                        }
                }else{
			if($peptideTypeHash{$exonstr} eq "Aluelement"){
        	                $uni{$_}{$pepCor}{$str} = '';
	                }else{
                        	$uni{$_}{$pepCor}{$str} = '' if $ele[9] !~ /exon/;
                	}
                }
	}
	# output
	foreach my $pepCor(keys %{$uni{$_}}){		
		my $len = scalar(keys %{$uni{$_}{$pepCor}});
		foreach my $s(keys %{$uni{$_}{$pepCor}}){
			say $pepCor,"\t",$s,"\t",$len;
		}
	}
}

system("rm -rf tmp/") if $rmtmp == 1;

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
				warn $line, "\n" if scalar(@arr) == 1;
				my ($chr, $strand) = @arr[0,4];
				my ($left, $jl) = (split /,/,$arr[5])[0,2];
				my ($jr, $right) = (split /,/,$arr[6])[0,2];
				my ($startTrim, $ORF) = @arr[7,8];
				$startTrim =~ s/startTrim://g;
				$ORF =~ s/ORF://g;
				$startTrim = $startTrim + $s;
				my $endPos = 1;
				# find comet peptide startPos
				my $crossSplice = "false";
				my ($a1, $a,  $a3) = ($left, 0, $jl);
			        my ($b1, $b2, $b3) = ($jr, 0, $right);
				my $startPos = $strand eq '+' ? ($a1+$startTrim*3+$ORF) : ($b3-$startTrim*3-$ORF);
				if($strand eq '+' && $startTrim*3 >= ($a3 - ($a1 + $ORF) + 1)){
			                $startPos = $b1 + $startTrim*3 - ($a3 - ($a1 + $ORF) + 1);			
			        }
			        if($strand eq '-' && $startTrim*3 >= ($b3 - $ORF - $b1 + 1)){
			                $startPos = $a3 - ($startTrim*3 - ($b3 - $ORF - $b1 + 1));
        			}
				my $endTag = $strand eq '+' ? ($startPos+length($pep)*3-1) : ($startPos-length($pep)*3+1);
				$crossSplice = 'true' if $strand eq '+' && $startPos <= $a3 && $endTag > $a3;
				$crossSplice = 'true' if $strand eq '-' && $startPos >= $b1 && $endTag < $b1;
				foreach my $ss((split/;/,$splice)){
					#warn $line,"\n" if $splice eq '17-52;0-17' && $pep =~ /YLLDLRNTSTPFKG/i;
					my ($splice_s, $splice_e) = split /-/, $ss;
					next if $e < $splice_s or $splice_e < $s;
					my $exonside = 'left';
					if(($splice_s =~ /^0/ && $seq =~ /^[a-z]/) or $splice_s !~ /^0/){
						$exonside = 'right';
					}
					if($crossSplice eq 'true'){ # peptide include splicing site
						if($strand eq '+'){
							$endPos = length($peptide) * 3 - ($jl - $startPos + 1) + $jr - 1;
							$endPos = length($peptide) * 3 + $startPos - 1 if $startPos > $jl;
							if($exonside eq 'right'){
								if($startPos > $jl){
									say OUT join "\t", ($chr,$startPos-1,$endPos,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
								}else{
									say OUT join "\t", ($chr,$jr-1,$endPos,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
								}
							}else{
								if($startPos > $jl){}else{
									say OUT join "\t", ($chr,$startPos-1,$jl,$id."_index:$s"."_left_".$pep."_".$_,0,$strand);
								}
							}
						}else{							
							$endPos = $jl + 1 - (length($peptide) * 3 - ($startPos - $jr + 1));
							$endPos = $startPos + 1 - length($peptide) * 3 if $startPos < $jr;
							if($exonside eq 'right'){
								if($startPos < $jr){
									say OUT join "\t", ($chr,$endPos-1,$startPos,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
								}else{
									say OUT join "\t", ($chr,$endPos-1,$jl,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
								}
							}else{
								if($startPos < $jr){}else{
									say OUT join "\t", ($chr,$jr-1,$startPos,$id."_index:$s"."_left_".$pep."_".$_,0,$strand);
								}
							}
						}
					}else{ # peptide include no splicing site
						if($strand eq '+'){
							$endPos = $startPos - 1 + length($peptide) * 3;
							if($exonside eq 'right'){
								say OUT join "\t", ($chr,$startPos-1,$endPos,$id."_index:$s"."_right_".$pep."_".$_,0,$strand);
							}else{
								say OUT join "\t", ($chr,$startPos-1,$endPos,$id."_index:$s"."_left_".$pep."_".$_,0,$strand);
							}
						}else{
							$endPos = $startPos + 1 - length($peptide) * 3;
							if($exonside eq 'right'){
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
my $str = "USAGE: $0 file.list(bin/tjScript result) bedfile(data/Alu.bed) SJ_outfile
Example:\n $0 result_v5/LCLs/HKgene/HK.nolocalFDR.0.05.list data/hk.sorted.bed 'data/SJ_out/*.SJ'\n";
return $str;
}
sub bed2hash {
	my $in = shift;
	die USAGE()."Please input bedfile\n" if !defined($in) or !-e $in;
	my %hash = ();
	open IN, $in;
	while(<IN>){
		chomp;
		my @arr = split;
		my $id = join ":",@arr[0,1,2,4,5];
		push @{$hash{$id}}, (join "_",@arr);
	}
	close IN;
	return \%hash;
}
sub SJ2hash {
	my $in = shift;
	
	my %hash = ();
	my %count = ();
	my @list = glob "$in";
	foreach my $f(@list){
		(my $basename = basename($f)) =~ s/\..*$|_(v|V).*$//g;
		open IN, $f;
		while(<IN>){
			chomp;
			my @arr = split;
			if(exists $hash{$basename}{$arr[1]}){
				$hash{$basename}{$arr[1]} = $hash{$basename}{$arr[1]} + $arr[6];
				$count{$basename}{$arr[1]} = $count{$basename}{$arr[1]} + 1;
			}else{
				$hash{$basename}{$arr[1]} = $arr[6];
				$count{$basename}{$arr[1]} = 1;
			}
			if(exists $hash{$basename}{$arr[2]}){
                                $hash{$basename}{$arr[2]} = $hash{$basename}{$arr[2]} + $arr[6];
				$count{$basename}{$arr[2]} = $count{$basename}{$arr[2]} + 1;
                        }else{
                                $hash{$basename}{$arr[2]} = $arr[6];
				$count{$basename}{$arr[2]} = 1;
                        }
			#$hash{$basename}{$arr[1]} = '';
			#$hash{$basename}{$arr[2]} = '';
		}
		close IN;
	}
	my %aboveOne = ();
	foreach my $sample(keys %hash){
		foreach my $site(keys %{$hash{$sample}}){
			if($hash{$sample}{$site}/$count{$sample}{$site} > 1){
				$aboveOne{$sample}{$site} = $hash{$sample}{$site};
			}
		}
	}
	return \%aboveOne;
	#return \%hash;
}
