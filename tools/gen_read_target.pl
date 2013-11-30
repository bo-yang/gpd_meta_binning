#!/usr/bin/perl

###
# gen_read_target.pl - This tool generate vectors of species for each read.
#
# Usage: gen_read_target.pl -i <input_reads> [-o <out_file>] 
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my @tmpFiles;

my $input;
my $output;
my $taxPath;
GetOptions('input|i=s'	 => \$input,
	   'output|o=s'  => \$output);

die("ERROR: no input file specied!\n") if (!$input);

my $basename=`basename $input`;
chomp($basename);
$output="$basename\_read_target.mat" if(!$output);

# Handle signals
my $CleanupDone   = 0;
$SIG{'ABRT'} = 'doCleanup';

# Ignore following signals.
$SIG{'HUP'}  = 'doCleanup';
$SIG{'INT'}  = 'doCleanup';
$SIG{'QUIT'} = 'doCleanup';

my $specs={}; # all species

open IN, "<$input" or die("ERROR: cannot open file $input!");

while (my $line = <IN>) {
	chomp $line;
	if($line =~ /^>/) {
		(my $read,my $loc,my $gi, my $gi_num, my $ref, my $ref_id, my $source_1)=split(/\|/,$line);
 	       $source_1 =~ s/ /_/g;
 	       $source_1 =~ s/,.*//g; # remove comma
 	       $source_1 =~ s/^_//g;  # leading _
 	       $source_1 =~ s/_chromosome.*//g;
 	       $source_1 =~ s/_plasmid.*//g;
 	       $source_1 =~ s/_incision.*//g;

 	       my $spec=$source_1;
       }
}


open OUT, ">$output" or die("ERROR: cannot create file $output!\n");

close OUT;

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

	close OUT;
	$CleanupDone = 1;
	exit 0;
} ## end sub doCleanup

sub Ignore_sig {
	# Do nothing
}
