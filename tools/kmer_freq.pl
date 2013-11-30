#!/usr/bin/perl

###
# kmer_filter.pl - This tool counts k-mer and separate reads based on k-mer distributions.
#
# Usage: kmer_filter.pl -i <input_genome> -o <out_file> -c <coverage> \
# 		-s <read_id_suffix> -l <read_length> -gl2 <gpd_lambda2>
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $kmerlen=16;
my $input;
my $output;

GetOptions('kmer|k=i'	=> \$kmerlen,
	   'input|i=s'  => \$input,
	   'output|o=s' => \$output
   );

die("ERROR: no input reads specied!\n") if (!$input);
my $basename=`basename $input`;
chomp($basename);
$output="$basename.$kmerlen-mer" if (!$output);

open IN, "<$input" or die("ERROR: cannot open file $input!");

my $kmerstat={};
my $read='';
my $kmersum=0;
while (my $line = <IN>) {
	chomp $line;
	if($line =~ /^>/) {
		if(length($read)!=0) {
			for(my $i=0;$i<length($read)-$kmerlen+1;++$i) {
				my $kmer=substr($read,$i,$kmerlen);
				if(not $kmerstat->{$kmer}) {
					$kmerstat->{$kmer}=1;
				} else {
					$kmerstat->{$kmer}++;
				}
				$kmersum++;
			}
			$read='';
		}
	} else {
		$read=$read.$line;
	} # end of if
}

if(length($read)!=0) {
	for(my $i=0;$i<length($read)-$kmerlen+1;++$i) {
		my $kmer=substr($read,$i,$kmerlen);
		if(not $kmerstat->{$kmer}) {
			$kmerstat->{$kmer}=1;
		} else {
			$kmerstat->{$kmer}++;
		}
		$kmersum++;
	}
	$read='';
}

# Record the k-mer statistics
open OUT, ">$output" or die("ERROR: cannot write file $output!");
#print OUT "k-mer\toccurrence\n";

foreach my $kmer (sort {$kmerstat->{$b} <=> $kmerstat->{$a}} keys %{$kmerstat}) {
#	my $ratio=$kmerstat->{$kmer}/$kmersum;
	print OUT "$kmer\t$kmerstat->{$kmer}\n";
}
