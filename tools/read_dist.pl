#!/usr/bin/perl

###
# kmer_dist.pl - This tool counts k-mer and calculate the Euclidean and KL distances.
#
# Usage: kmer_dist.pl -k <kmer> -i <input_genome> -o <out_file> -c <coverage> 
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $kmerlen=16;
my @input;
my $output;

GetOptions('kmer|k=i'	=> \$kmerlen,
	   'input|i=s'  => \@input,
	   'output|o=s' => \$output
   );

@input = split(/,/,join(',',@input));

die("ERROR: no input reads specied!\n") if (!$input);
my $basename=`basename $input`;
$output="$kmerlen-mer_$basename" if (!$output);

open IN, "<$input" or die("ERROR: cannot open file $input!");

my $kmerstat={};
my $read='';
my $read_id='';
#my $kmersum=0;
while (my $line = <IN>) {
	chomp $line;
	if($line =~ /^>/) {
		if(length($read)!=0 && length($read_id)!=0) {
			$kmerstat->{$read_id}->$kmersum=0;
			for(my $i=0;$i<length($read)-$kmerlen+1;++$i) {
				my $kmer=substr($read,$i,$kmerlen);
				if(not $kmerstat->{$read_id}->{$kmer}) {
					$kmerstat->{$read_id}->{$kmer}=1;
				} else {
					$kmerstat->{$read_id}->{$kmer}++;
				}
				$kmerstat->{$read_id}->$kmersum++;
			}
			$read='';
		}
		($read_id,my $sources,my $errors,my $source_1)=split(/\|/,$line);
	} else {
		$read=$read.$line;
	} # end of if
}
# Don't forget the last read!
if(length($read)!=0 && length($read_id)!=0) {
	$kmerstat->{$read_id}->$kmersum=0;
	for(my $i=0;$i<length($read)-$kmerlen+1;++$i) {
		my $kmer=substr($read,$i,$kmerlen);
		if(not $kmerstat->{$read_id}->{$kmer}) {
			$kmerstat->{$read_id}->{$kmer}=1;
		} else {
			$kmerstat->{$read_id}->{$kmer}++;
		}
		$kmerstat->{$read_id}->$kmersum++;
	}
	$read='';
}

# Factorize each k-mer frequency
foreach my $id (keys %$kmerstat) {
	$kmerstat->{$read_id}->{$kmer} /= $kmerstat->{$read_id}->$kmersum;
}

# Calculate the 



# Record the k-mer statistics
open OUT, ">$output" or die("ERROR: cannot write file $output!");
#print OUT "k-mer\toccurrence\n";

foreach my $kmer (sort {$kmerstat->{$b} <=> $kmerstat->{$a}} keys %{$kmerstat}) {
	my $ratio=$kmerstat->{$kmer}/$kmersum;
	print OUT "$kmer\t$ratio\n";
}
