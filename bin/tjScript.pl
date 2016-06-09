#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

#tjRNAJunction($ARGV[0]);
tjPercolaotr($ARGV[0]);
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
	my $in = shift;

	(my $dir = 'result_v5/LCLs/Alu_all/') =~ s/\/$//g;
        my $fdr = 0.05;
        my $bedfile = '/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/Ensembl_Alu_25bp_0.5.unique.sorted.bed';
	my @rawfiles = `ls $dir/cometout/`;

        my %hash = ();
        my $n = 1;
        foreach(@rawfiles){
                chomp;
                $hash{$n} = $_;
                $n++;
        }
	die $in,"\n" if !exists $hash{$in};

	my $junction = $dir."/junctionPep/$hash{$in}.fa";
	(my $percolator = $dir."/cometout/$hash{$in}/$fdr") =~ s/\/$//g;
	#say "/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/pep_percolator.pl $junction $percolator $bedfile";
	system("/u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/bin/pep_percolator_fast.pl $junction $percolator $bedfile");
	#system("rm -rf tmp");
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
GM18855.rna
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
GM19143.rna
GM19144.rna
GM19147.rna
GM19152.rna
GM19153.rna
GM19160.rna
GM19172.rna
GM19192.rna
GM19193.rna
GM19200.rna
GM19204.rna
GM19207.rna
GM19209.rna
GM19222.rna
GM19257.rna

