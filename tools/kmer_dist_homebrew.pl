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

GetOptions('kmer|k=i'	=> \$kmerlen,
	   'input|i=s'  => \$input
   );

die("ERROR: no input reads specied!\n") if (!$input);

open IN, "<$input" or die("ERROR: cannot open file $input!");

my $kmerstat={};
while (my $line = <IN>) {
	chomp $line;
	if($line =~ /^>/) {
		# Record read ID?
	} else {
		for(my $i=0;$i<length($line)-$kmerlen+1;++$i) {
			my $kmer=substr($line,$i,$kmerlen);
			if(not $kmerstat->{$kmer}) {
				$kmerstat->{$kmer}=1;
			} else {
				$kmerstat->{$kmer}++;
			}
		}
	} # end of if
}

# Record the k-mer statistics
open OUT, ">$input.kmer$kmerlen.out" or die("ERROR: cannot write file $input.out!");
print OUT "k-mer\toccurrence\n";

foreach my $kmer (sort {$kmerstat->{$b} <=> $kmerstat->{$a}} keys %{$kmerstat}) {
	print OUT "$kmer\t$kmerstat->{$kmer}\n";
}
