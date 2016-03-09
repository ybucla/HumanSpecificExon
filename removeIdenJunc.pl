#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

my $in = shift or die "$0 fasta";

die "No such file found!\n" if !-e $in;

my %seq = ();
my $head = '';
open IN, $in;
while(<IN>){
	chomp;
	if(/^>/){
		($head = $_) =~ s/^>//g;
		$seq{$head} = '';
	}else{
		$seq{$head} = $seq{$head}.$_; 
	}
}
close IN;

# do blast search
my @blastout = `blastall -p blastp -i $in -d data/blastDb/uniprot_1ident_58896.fasta -m8 -b1`;
my %hash = ();
foreach(@blastout){
	chomp;
	my ($id,$uid, $identity,$alignlen) = (split /\s+/)[0,1,2,3];
	next if !exists $seq{$id};
	my $orilen = length($seq{$id});
	if($orilen == $alignlen && $identity >= 100){		
		delete($seq{$id});
		$hash{$id} = $uid;
	}
}

# output merged database for ms search
my $out = "massDb/merge_".basename($in);
my @proteom = `cat data/blastDb/uniprot_1ident_58896.fasta`;
open OUT, ">$out";
print OUT $_ for @proteom;
print OUT ">",$_,"\n",$seq{$_},"\n" for keys %seq;
close OUT;

mkdir "massDb/identity" if !-e "massDb/identity";
open OUTB, ">massDb/identity/identify_".basename($in);
print OUTB $_,"\t",$hash{$_},"\n" for keys %hash;
close OUTB;



__DATA__
GM18486.rna
GM18498.rna
GM18499.rna
GM18501.rna
GM18502.rna
GM18504.rna
GM18505.rna
GM18507.rna
GM18508.rna
GM18510.rna
GM18511.rna
GM18516.rna
GM18517.rna
GM18519.rna
GM18520.rna
GM18522.rna
GM18523.rna
GM18852.rna
GM18853.rna
GM18855.rna
GM18856.rna
GM18858.rna
GM18861.rna
GM18862.rna
GM18870.rna
GM18909.rna
GM18912.rna
GM18913.rna
GM18916.rna
GM19093.rna
GM19098.rna
GM19099.rna
GM19101.rna
GM19102.rna
GM19108.rna
GM19114.rna
GM19116.rna
GM19119.rna
GM19127.rna
GM19128.rna
GM19130.rna
GM19131.rna
GM19137.rna
GM19138.rna
GM19140.rna
GM19141.rna
GM19143.rna
GM19144.rna
GM19147.rna
GM19152.rna
GM19153.rna
GM19159.rna
GM19160.rna
GM19171.rna
GM19172.rna
GM19190.rna
GM19192.rna
GM19193.rna
GM19200.rna
GM19201.rna
GM19204.rna
GM19207.rna
GM19209.rna
GM19210.rna
GM19222.rna
GM19225.rna
GM19238.rna
GM19239.rna
GM19257.rna

