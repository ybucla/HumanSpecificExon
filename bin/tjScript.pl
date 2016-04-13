#! /u/home/y/ybwang/perl

# get floss value for each region, against with coding region of known gene

use strict;
use warnings;
use File::Basename;
use Data::Table;
use 5.010;

#tjRNAJunction($ARGV[0]);
tjPercolaotr($ARGV[0]);

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
	my $in = shift;

	my $dir = 'result_v5/HSE/';
	my $fdr = 0.01;
	my $i = 0;
        my @rawfiles = `ls $dir/cometout/`;
        foreach(@rawfiles){
                chomp;
		next if /PJ_/;
		$i += 1;
		next if defined($in) &&  $i != $in;
		
	        my $percolatorfile = $dir.'/cometout/'.$_.'/'.$fdr.'/*.'.$fdr.'.pep.percolator';

        	my $reg = 'chr';
		#my $identityfile = $dir.'/massDb/identity/identify_'.$_.'.fa';
	        #my @identityhead = `cut -f 2 $identityfile | sort | uniq`;
        	#foreach(@identityhead){
	        #        chomp;
                #	$reg = $reg.'|'.$_;
        	#}

        	my @n1 = `awk -F '\\t' '\$3 < $fdr' $percolatorfile | cut -f 6 | sort | uniq`;
	        my @n2 = `awk -F '\\t' '\$3 < $fdr' $percolatorfile | awk -F '\\t' '\$6 ~ /$reg/' | cut -f 6 | sort | uniq`;
        	my @m1 = `awk -F '\\t' '\$3 < $fdr' $percolatorfile | sort | uniq`;
	        my @m2 = `awk -F '\\t' '\$3 < $fdr' $percolatorfile | awk -F '\\t' '\$6 ~ /$reg/' | sort | uniq`;
	        print "======", $_, "\n";
	        print "Protein number:\t", scalar(@n1), "\n";
        	#print "Junction number:\t", $n2[0];
	        #print "Pep number:\t", $m1[0];
        	print "Junction Pep number:\t", scalar(@n2),"\n";
		
		print $_ for @n2;
        }
}
