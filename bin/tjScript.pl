#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

#tjRNAJunction($ARGV[0]);
tjPercolaotr();
#tjPercolaotr_TISSUE($ARGV[0]);

sub tjRNAJunction {
	my $in = shift;
	
	my $n = 1;
	my @rawfiles = `ls cometout/`;
        foreach(@rawfiles){
                chomp;
		next if defined($in) && $n != $in;
		my $bamfile = "RNA/".$_."/Aligned.out.sorted.bam";
		$bamfile = "/u/home/y/ybwang/scratch/Mathias_Wilhelm/rna/".$_."/star_out/Aligned.out.sorted.bam" if $_ =~ /lymph/;
		say $_;
		my @result = `samtools view $bamfile | awk -F '\\t' '\$6 ~ /N/' | cut -f 1 | sort | uniq | wc -l`;
                print $result[0];
		
		$n++;
        }
}

sub tjPercolaotr {
	#my $in = shift;

	(my $dir = 'result_v5/LCLs/HKgene/') =~ s/\/$//g;;
	my $fdr = 0.05;
	my $bedfile = '/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/hk.sorted.bed';

	my $i = 0;
        my @rawfiles = `ls $dir/cometout/`;
        foreach(@rawfiles){
                chomp;
		next if /PJ_/;
		my $junction = $dir."/junctionPep/$_.fa";
		(my $percolator = $dir."/cometout/$_/$fdr") =~ s/\/$//g;
		#say "/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/pep_percolator.pl $junction $percolator $bedfile";
		system("/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/pep_percolator.pl $junction $percolator $bedfile");
        }
	system("rm -rf tmp");
}

sub tjPercolaotr_TISSUE {
        my $in = shift; # lymph
	die "USAGE:\t$0 tissuename(i.e. liver)\n" if !defined($in);

        my $dir = 'result_v5/TISSUE/Alu/cometout/'.$in;
        my $fdr = 0.05;
        my $bedfile = '/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/Ensembl_Alu_25bp_0.5_CDS_overlap.unique.sorted.bed';
	(my $junction = $dir) =~ s/cometout/junctionPep/;
	$junction .= '.fa';
	my $percolator = $dir.'/'.$fdr;
	
	#say "/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/pep_percolator.pl $junction $percolator $bedfile";
	system("/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/pep_percolator.pl $junction $percolator $bedfile");
	#system("rm -rf tmp");
}
