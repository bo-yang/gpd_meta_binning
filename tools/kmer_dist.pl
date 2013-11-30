#!/usr/bin/perl

###
# kmer_dist.pl - This tool counts k-mer and separate reads based on k-mer distributions.
#
# Usage: kmer_dist.pl -k <kmer> -i <genome1,genome2> 
# 		
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $kmerlen=3;
my @input;
my $output;

GetOptions('kmer|k=i'	=> \$kmerlen,
	   'input|i=s'  => \@input,
	   'output|o=s' => \$output
   );

@input = split(/,/,join(',',@input));

die("ERROR: no input reads specied!\n") if (!@input);

my $kmerstat={};
my $speciesfile=''; # record the whole species file number
my $maxkmer=4**$kmerlen;
foreach my $file (@input) {
	open IN, "<$file" or die("ERROR: cannot open file $file!");
	
	my $read='';
	my $kmersum=0;
	while (my $line = <IN>) {
		chomp $line;
		if($line =~ /^>/) {
			if(length($read)!=0) {
				for(my $i=0;$i<length($read)-$kmerlen+1;++$i) {
					my $kmer=substr($read,$i,$kmerlen);
					if(not $kmerstat->{$file}->{$kmer}) {
						$kmerstat->{$file}->{$kmer}=1;
					} else {
						$kmerstat->{$file}->{$kmer}++;
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
			if(not $kmerstat->{$file}->{$kmer}) {
				$kmerstat->{$file}->{$kmer}=1;
			} else {
				$kmerstat->{$file}->{$kmer}++;
			}
			$kmersum++;
		}
		$read='';
	}
	
	my $basename=`basename $file`;
	$output="$kmerlen-mer_$basename";
	open OUT, ">$output" or die("ERROR: cannot write file $output!");

	my $kmernum=0;
	# Facterize the k-mer statistics
	foreach my $kmer (keys %{$kmerstat->{$file}}) {
		$kmernum++;
		$kmerstat->{$file}->{$kmer} /= $kmersum;
		print OUT "$kmer\t$kmerstat->{$file}->{$kmer}\n";
	}

	if($kmernum>$maxkmer){	# Detect the file with max number of k-mers
		$speciesfile=$file;
		$maxkmer=$kmernum;
	}
#	print "$file: $kmernum kmers\n"; # TEST ONLY

	close OUT;
	close IN;
}

#print "Species file: $speciesfile\n"; # TEST ONLY

my $eucsum=0.0; # Euclidean distance
my $klsum=0.0;  # KL distance
foreach my $kmer (keys %{$kmerstat->{$speciesfile}}) {
	if(not $kmerstat->{$input[0]}->{$kmer}) {
		$kmerstat->{$input[0]}->{$kmer}=0.0000000001;
	}
	if(not $kmerstat->{$input[1]}->{$kmer}) {
		$kmerstat->{$input[1]}->{$kmer}=0.0000000001;
	}

	$eucsum += ($kmerstat->{$input[0]}->{$kmer} - $kmerstat->{$input[1]}->{$kmer})**2;
	$klsum += $kmerstat->{$input[0]}->{$kmer}*log($kmerstat->{$input[0]}->{$kmer}/$kmerstat->{$input[1]}->{$kmer});
}
$eucsum = sqrt($eucsum);

print "Euclidean distance: $eucsum\n";
print "KL distance: $klsum\n";
