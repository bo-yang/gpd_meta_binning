#!/usr/bin/perl

###
# bin_reads.pl - This tool clustering reads based on clustered k-mers.
#
# Usage: bin_reads.pl -i <input_reads> -c <clustered_k-mers> -o <out_file>
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $input;
my $cluster;
my $output;

GetOptions('cluster|c=s'=> \$cluster,
	   'input|i=s'  => \$input,
	   'output|o=s' => \$output
   );

die("ERROR: no input reads specied!\n") if (!$input);
my $basename=`basename $input`;
chomp($basename);
$output="$basename.out" if (!$output);

# Read k-mer clusters
open CLIN, "<$cluster" or die("ERROR: cannot open file $cluster!");
my $kmerclst={};
my $numclst={};	# number of different clusters
my $kmerlen=0;
while (my $line = <CLIN>) {
	chomp $line;
	(my $kmer,my $clst)=split("\t",$line);
	$kmerclst->{$kmer}=$clst;
	if(not $numclst->{$clst}){
		$numclst->{$clst}=1;
	} else {
		$numclst->{$clst}++;
	}
	$kmerlen=length($kmer);
}
close CLIN;

# Classify raw reads based on clustered k-mers.
open IN, "<$input" or die("ERROR: cannot open file $input!");

my $header='';
my $body='';
my $kmercnt={};

while (my $line = <IN>) {
	chomp $line;

	if($line =~ /^>/) {
		if(length($body)!=0) {
			# count number of unique k-mers
			for(my $i=0;$i<length($body)-$kmerlen+1;++$i) {
				my $kmer=substr($body,$i,$kmerlen);
				if(not $kmercnt->{$kmer}) {
					$kmercnt->{$kmer}=1;
				} else {
					$kmercnt->{$kmer}++;
				}
			} # end of for

			# get the k-mer of max occurrence
			my @sortedkmer=sort {$kmercnt->{$b} <=> $kmercnt->{$a}} keys %{$kmercnt};
			my $dormkmer=$sortedkmer[0];
			my $rdclst=$kmerclst->{$dormkmer};
			#print "Best k-mer: $dormkmer => $rdclst(occur: $kmercnt->{$dormkmer})\n"; # TEST ONLY

			# write this read to file
			open OUT, ">>$output.$rdclst" or die("ERROR: cannot write file $$output.$rdclst!");
			print OUT "$header\n";
			print OUT "$body\n";
			close OUT;

			$body='';
			$kmercnt={};
		}

		$header=$line; # Remember header info

	} else {
		$body.=$line;
	} # end of if
}

if(length($body)!=0){
	# count number of unique k-mers
	for(my $i=0;$i<length($body)-$kmerlen+1;++$i) {
		my $kmer=substr($body,$i,$kmerlen);
		if(not $kmercnt->{$kmer}) {
			$kmercnt->{$kmer}=1;
		} else {
			$kmercnt->{$kmer}++;
		}
	} # end of for

	# get the k-mer of max occurrence
	my @sortedkmer=sort {$kmercnt->{$b} <=> $kmercnt->{$a}} keys %{$kmercnt};
	my $dormkmer=$sortedkmer[0];
	my $rdclst=$kmerclst->{$dormkmer};

	# write this read to file
	open OUT, ">>$output.$rdclst" or die("ERROR: cannot write file $$output.$rdclst!");
	print OUT "$header\n";
	print OUT "$body\n";
	close OUT;

	$body='';
	$kmercnt={};
}

close IN;

