#!/usr/bin/perl

###
# kmeans.pl - This tool clustering metagenomic reads by K-Means method.
#
# Usage: kmeans.pl -i <input_genome> -o <out_file> -k <kmer_len>
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';
use Algorithm::KMeans;

my $kmerlen=4;
my $input;
my $output;

GetOptions('input|i=s'  => \$input,
	   'output|o=s' => \$output,
	   'kmer|k=i'   => \$kmerlen
   );

die("ERROR: no input reads specied!\n") if (!$input);
die("ERROR: unsupported k-mer length!\n") if ($kmerlen>5 || $kmerlen<3);

my $basename=`basename $input`;
chomp($basename);
$output="$basename.cluster" if (!$output);

# Generate mask: "N1111...1"
my $mask='N';
for(my $i=0;$i<4**$kmerlen;++$i) {
	$mask=$mask.'1';
}

# Construct an instance of the clusterer.
my $clusterer=Algorithm::KMeans->new(
	datafile	=> $input,
	mask		=> $mask,
	K		=> 3,
	cluster_seeding	=> 'smart',
	terminal_output	=> 1,
	debug		=> 0,
	write_clusters_to_files	=> 1,
);

$clusterer->read_data_from_file();
$clusterer->kmeans();

