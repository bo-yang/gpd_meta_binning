#!/usr/bin/perl

###
# gen_kmer_mat.pl - This tool generates all possible k-mers, counts k-mer per read 
# 			and write these k-mer-occur-per-read into file.
#
# Usage: gen_kmer_mat.pl -i <input_genome> -o <out_file> 
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

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
$output="$basename.$kmerlen-mer.mat" if (!$output);

open IN, "<$input" or die("ERROR: cannot open file $input!");

my $kmerstat={};

# generate all possible k-mers
my @nt=('A','C','G','T');

if($kmerlen==3) {
	for(my $i=0;$i<4;$i++) {
		for(my $j=0;$j<4;$j++) {
			for(my $k=0;$k<4;$k++) {
				$kmerstat->{"$nt[$i]$nt[$j]$nt[$k]"}=1;
			}
		}
	} # end of for
} elsif($kmerlen==4) {
	for(my $i=0;$i<4;$i++) {
		for(my $j=0;$j<4;$j++) {
			for(my $k=0;$k<4;$k++) {
				for(my $l=0;$l<4;$l++) {
					$kmerstat->{"$nt[$i]$nt[$j]$nt[$k]$nt[$l]"}=1;
				}
			}
		}
	} # end of for
} elsif($kmerlen==5) {
	for(my $i=0;$i<4;$i++) {
		for(my $j=0;$j<4;$j++) {
			for(my $k=0;$k<4;$k++) {
				for(my $l=0;$l<4;$l++) {
					for(my $m=0;$m<4;$m++) {
						$kmerstat->{"$nt[$i]$nt[$j]$nt[$k]$nt[$l]$nt[$m]"}=1;
					}
				}
			}
		}
	} # end of for
}

# Record the k-mer statistics
open OUT, ">$output" or die("ERROR: cannot write file $output!");
my @allkmers=keys %$kmerstat;
#print OUT "@allkmers\n";

# Count k-mer statistics
my $read='';
my $readID;
while (my $line = <IN>) {
	chomp $line;
	if($line =~ /^>/) {
		(my $ID,my $junk)=split(/\|/,$line);
		$ID =~ s/^>r//g;
		$ID =~ s/\.//g;
		$ID =~ s/ //g;

		if(length($read)!=0) {
			# count k-mers for this read
			for(my $i=0;$i<length($read)-$kmerlen+1;++$i) {
				my $kmer=substr($read,$i,$kmerlen);
				if(not $kmerstat->{$kmer}) {
					# noise, just ignore it
					print "Unknown k-mer: $kmer\n";
				} else {
					$kmerstat->{$kmer}++;
				}
			}

			# write k-mers and reset all variables
			print OUT "$readID\t"; # symbolic ID.
			$read='';
			my $cnt='';
			foreach my $kmer (keys %$kmerstat) {
				$kmerstat->{$kmer}--; # minus the default 1 occurrence
				$cnt=$cnt."$kmerstat->{$kmer} "; 
				$kmerstat->{$kmer}=1;
			}
			print OUT "$cnt\n";
		}
		$readID=$ID;
	} else {
		$read=$read.$line;
	} # end of if
}

# Don't forget the last read
if(length($read)!=0) {
	# count k-mers for this read
	for(my $i=0;$i<length($read)-$kmerlen+1;++$i) {
		my $kmer=substr($read,$i,$kmerlen);
		if(not $kmerstat->{$kmer}) {
			# noise, just ignore it
			print "Unknown k-mer: $kmer\n";
		} else {
			$kmerstat->{$kmer}++;
		}
	}

	# write k-mers and reset all variables
	print OUT "$readID\t"; # symbolic ID.
	$read='';
	my $cnt='';
	foreach my $kmer (keys %$kmerstat) {
		$kmerstat->{$kmer}--; # minus the default 1 occurrence
		$cnt=$cnt."$kmerstat->{$kmer} ";
		$kmerstat->{$kmer}=1;
	}
	print OUT "$cnt\n";
}


close OUT;
close IN;
