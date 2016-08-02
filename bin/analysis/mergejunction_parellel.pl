#! /u/home/y/ybwang/perl
# merge junction sequence for single individuals with cdhit (used for human tissue proteome)
use strict;
use warnings;
use File::Basename;
use Data::Table;
use Thread;
use 5.010;

my $in = shift or die "ERROR: $0 liver\n"; # 'liver'
my $fasta = "single_junctionPep/$in*";

# read fasta file
my $seq = readSeq($fasta);

# read list file and write each junctino seq
my %junction;
my $head = '';
foreach my $head(keys %{$seq}){
	my @arr = split /,|_/,$head;
        my $id = join "_", @arr[0,1,2,4];
        $junction{$id}{$head} = $seq->{$head};
}


my $tmpdir = 'tmp_'.$in;
mkdir $tmpdir if !-e $tmpdir;

# thread start
my %result_hash :shared = ();
my $index :shared = 0;

my @id = keys %junction;
my @threads;
my $tempcount = 0;
foreach(1..1){
	$threads[$tempcount]=Thread->new(\&run,\@id,$_,$tmpdir,\%junction);
	$tempcount++;
}
foreach my $thread (@threads) {
	my $r = $thread->join();
}

open ROUT, ">$in.fa";
say ROUT ">",$_,"\n",$result_hash{$_} for keys %result_hash;
close ROUT;

system("rm -rf $tmpdir/");
warn "# finished!\n";

# --sub--
sub readSeq {
	my $f = shift;
	my @list = glob $f;
	my %seq = ();
	foreach(@list){
		(my $basename = basename($_)) =~ s/\..*$//g;
		open IN, $_;
		my $id = '';
		while(<IN>){
			chomp;
			if(/>/){
				($id = $_) =~ s/>//g;
			}else{
				#$seq{$id."_".$basename} = $_;
				$seq{$id} = $_;
			}
		}
		close IN;
	}
	return \%seq;
}

sub cdhit {
	my $in = shift;
	my $out = $in.'.cdhit';
	my $bin = '/u/home/y/ybwang/nobackup-yxing/program/cd-hit-v4.6.5-2016-0304/cd-hit';
	#say "$bin -i $in -d 0 -o $out -c 1.0 -n 5 -l 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0";
	my @r = `$bin -i $in -d 0 -o $out -c 1.0 -n 5 -l 5 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0`;	
}

sub run {
	my $rarr = shift;
	my $n = shift;
	my $tmpdir = shift;
	my $junction = shift;

	my @arr = @{$rarr};
	while( (my $i = getindex()) < scalar(@arr)){	
		warn "Thread-$n: $i\n";
		my $k1 = $arr[$i];
		if($k1 =~ /chr10_100011459_100012110_/){
			#foreach my $head(keys %{$junction->{$k1}}){
			#	say $head,"\t",$junction->{$k1}{$head};
			#}
		}
		my $arr = cluster($junction->{$k1});
		{
			lock(%result_hash);
			foreach(@{$arr}){
				$result_hash{$_} = $seq->{$_};
			}
		}
		#my $k1 = $arr[$i];
		#open OUT, ">$tmpdir/".$k1;
		#foreach my $k2(keys %{$junction->{$k1}}){
		#	say OUT ">",$k2,"\n",$junction->{$k1}{$k2};
		#}
		#close OUT;
		## cdhit search
		#cdhit("$tmpdir/".$k1);
		#my @result = `cat $tmpdir/$k1.cdhit.clstr`;
		#{
		#	lock(%result_hash);
		#	foreach(@result){
        	#		chomp;
		#	        if(/>(.*)\.\.\.\s+\*/){				
                #			$result_hash{$1} = $seq->{$1};
        	#		}
		#	}
		#}
		#system("rm $tmpdir/$k1*");
	}
}

sub cluster {
	my $hash = shift;
	my @key = keys %{$hash};
	my %redandunt = ();
	foreach my $i(0..$#key){
		next if exists $redandunt{$i};
		my $seqi = uc($hash->{$key[$i]});
		foreach my $j($i+1..$#key){
			next if exists $redandunt{$j};
			my $seqj = uc($hash->{$key[$j]});
			if($seqi eq $seqj){
				$redandunt{$j} = '';
			}else{
				$redandunt{$j} = '' if index($seqi, $seqj) != -1;
				$redandunt{$i} = '' if index($seqj, $seqi) != -1;
			}
		}
	}
	my @ids = ();
	foreach(0..$#key){
		push @ids, $key[$_] if !exists $redandunt{$_};
	}
	return \@ids;
}

sub getindex {

	lock($index);

	return $index++;
}
