#!/usr/bin/perl

###
# gpdbl.pl - This tool calls GPD program and BLASTN to do metagenomic binning.
#
# Usage: gpdbl.pl -i <input_reads> [-o <out_file>] [-db <blast_db>] [-t <taxonomy_db>]
#	input_reads: reads of genome to be queried, in FASTA format.
#	output_file: the name of output files
#	blast_db:    BLAST databse to be used
#	taxonomy_db: path to taxonomy DBs
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
	   'output|o=s'  => \$output
   	   );

die("ERROR: no input reads specied!\n") if (!$input);

# Handle signals
my $CleanupDone   = 0;
$SIG{'ABRT'} = 'doCleanup';

# Ignore following signals.
$SIG{'HUP'}  = 'doCleanup';
$SIG{'INT'}  = 'doCleanup';
$SIG{'QUIT'} = 'doCleanup';

my $genMycoplasma=580076; # length of genome
my $genBuchnera=615980; # length of genome

###########################################################################################
# 
# Count species and get the p-values
#
###########################################################################################
my $read_len={};
$read_len->{'Mycoplasma'}=0;
$read_len->{'Buchnera'}=0;

my $isSpec={};
$isSpec->{'Mycoplasma'}=0;
$isSpec->{'Buchnera'}=0;

open IN, "<$input" or die("ERROR: cannot open file $input!");

while (my $line = <IN>) {
	chomp $line;
	if($line =~ /^>/) {
		(my $read,my $sources,my $errors,my $source_1)=split(/\|/,$line);
		if($source_1 =~ /\"(.*)\"/g) {
			my $spec=$1;

			if($spec =~ /Mycoplasma/) {
				$isSpec->{'Mycoplasma'}=1;
				$isSpec->{'Buchnera'}=0;
			} 
			if($spec =~ /Buchnera/) {
				$isSpec->{'Mycoplasma'}=0;
				$isSpec->{'Buchnera'}=1;
			}
		}
	} else {
		if($isSpec->{'Mycoplasma'}) {
			$read_len->{'Mycoplasma'}+=length($line);
		} 
		if($isSpec->{'Buchnera'}) {
			$read_len->{'Buchnera'}+=length($line);
		}
		$isSpec->{'Mycoplasma'}=0;
		$isSpec->{'Buchnera'}=0;
	}
}

close IN;

print "Total length of Mycoplasma: $read_len->{'Mycoplasma'}.\n";
print "Total length of Buchnera: $read_len->{'Buchnera'}.\n";
my $coverage=$read_len->{'Mycoplasma'}/$genMycoplasma;
print "Coverage(Mycoplasma)=$coverage.\n";
my $coverage=$read_len->{'Buchnera'}/$genBuchnera;
print "Coverage(Buchnera)=$coverage.\n";

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
