#!/usr/bin/perl

###
# metabin.pl - This tool calls GPD program,k-means clustering and GC content to do metagenomic binning.
#
# Usage: metabin.pl -i <input_reads> [-o <out_file>] [-t <taxonomy_db>]
#	input_reads: reads of genome to be queried, in FASTA format.
#	output_file: the name of output files
#	taxonomy_db: path to taxonomy DBs
###

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my @tmpFiles;

my $input;
my $output;
my $blastDB;
my $taxPath;
my $kmerLen=-1;
GetOptions('input|i=s'	 => \$input,
	   'output|o=s'  => \$output,
   	   'taxonomy|t=s'=> \$taxPath,
   	   'kmer|k=i'	 => \$kmerLen);

die("ERROR: no input reads specied!\n") if (!$input);

$input=abs_path($input);
my $input_file=`basename $input`;
chomp $input_file;

if($kmerLen != -1) {
	if(!$output) {
		$output="$input_file"."_kmer"."$kmerLen".".out" ;
	} else {
		$output="$output"."_kmer"."$kmerLen";
	}
}
$output="$input_file.out" if(!$output);

my $logFile="$Bin/$output\_log";
open LOG, ">$logFile" or die("Can't open $logFile for writing.\n");

# Handle signals
my $CleanupDone   = 0;
$SIG{'ABRT'} = 'doCleanup';

# Ignore following signals.
$SIG{'HUP'}  = 'doCleanup';
$SIG{'INT'}  = 'doCleanup';
$SIG{'QUIT'} = 'doCleanup';

###########################################################################################
#
# Call GPD program to classify genome reads
#
###########################################################################################
print "\nBinning reads by Generalised Poisson Distribution...\n";
print LOG "\nBinning reads by Generalised Poisson Distribution...\n";

die "ERROR: failed to find executable gpd under $Bin!\n" if !-x "$Bin/gpd";

my $gpdout="$Bin/.$output";
$ENV{'LD_LIBRARY_PATH'}="$ENV{'LD_LIBRARY_PATH'}:/usr/local/lib";
my $gpdcmd;
if($kmerLen != -1) {
	$gpdcmd="$Bin/gpd -input $input -output $gpdout -kmer_len $kmerLen -auto";
} else {
	$gpdcmd="$Bin/gpd -input $input -output $gpdout -auto"
}
system("$gpdcmd"); # run GPD
my $bestNumofBins=get_bin_no("$gpdout.score");
print "Best number of bins: $bestNumofBins\n";
print LOG "$bestNumofBins species have been identified in GPD.\n";

my $realNumofBins=$bestNumofBins;
my @gpdOutputs;
push(@tmpFiles,"$gpdout.score");
for(my $cnt=0;$cnt<$bestNumofBins;++$cnt) {
	if(-e "$gpdout.$cnt") {
		push(@tmpFiles,"$gpdout.$cnt"); # record temp files
		push(@gpdOutputs,"$gpdout.$cnt"); 
	} else {
		$realNumofBins--;
	}
}

my $date = `date`;
chomp $date;
print "done.  [$date]\n\n";

###########################################################################################
#
# Call K-means to further classify genome reads
#
###########################################################################################


###########################################################################################
#
# Check the GC content distribution
#
###########################################################################################



#####################################################
#
#	SUBROUTINES
#
#####################################################
sub get_bin_no {
	my $score_file=shift;
	my $bin_no;

	# TODO: calculate scores of GPD classification
	open IN, "<$score_file" or die("Can't open $score_file for reading.\n");
	while ( my $line = <IN> ) {
		chomp $line;
		if($line =~ /Number|number/) { 
			my @fields = split(/:/, $line);
			my $lable = $fields[0];
			$bin_no = $fields[1];
			last;
		}
	}

	close IN;
	return $bin_no;	
}

#########################################################
#               doCleanup
#########################################################
sub doCleanup {
	if($CleanupDone == 1) {
		exit;
	}

	# Ignore the following signals.
	$SIG{'HUP'}  = 'Ignore_sig';
	$SIG{'INT'}  = 'Ignore_sig';
	$SIG{'QUIT'} = 'Ignore_sig';
	$SIG{'ILL'}  = 'Ignore_sig';
	$SIG{'PIPE'} = 'Ignore_sig';
	$SIG{'ALRM'} = 'Ignore_sig';
	$SIG{'USR1'} = 'Ignore_sig';
	$SIG{'USR2'} = 'Ignore_sig';
	$SIG{'STOP'} = 'Ignore_sig';
	$SIG{'TSTP'} = 'Ignore_sig';
	$SIG{'CONT'} = 'Ignore_sig';
	$SIG{'TTIN'} = 'Ignore_sig';
	$SIG{'TTOU'} = 'Ignore_sig';
	$SIG{'LOST'} = 'Ignore_sig';
	$SIG{'CHLD'} = 'Ignore_sig';

	# delete temp files
	foreach my $file (@tmpFiles) {
		system("rm -f $file");
	}

	close LOG;
	$CleanupDone = 1;
	exit 0;
} ## end sub doCleanup

sub Ignore_sig {
	# Do nothing
}
