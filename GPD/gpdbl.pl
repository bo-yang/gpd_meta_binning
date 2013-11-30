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
	   'db|b=s'	 => \$blastDB,
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
# Read organism names so we can generate results files comparable to what BLAST will output.
# 
###########################################################################################

print "\nLoading taxonomy...\n";
print LOG "\nLoading taxonomy...\n";

# Accession scan 1: RefSeq orgs.
my $accFile = "$taxPath/.0_accessionMap/accessionMap.txt";
open IN, "<$accFile" or die("Can't open $accFile for reading.\n");
my $speciesDirName = {};

while ( my $line = <IN> ) {
   chomp $line;
   (my $orgName, my $prefix, my $seqType, my $desc) = split(/\t/, $line);
   $speciesDirName->{$prefix} = $orgName;
}

close IN;

###########################################################################################
# 
# Load full taxonomic metadata for all organisms in the database.
# 
###########################################################################################

# Taxonomy scan 1: RefSeq organisms.
my $taxFile = "$taxPath/.3_parsedTaxData/distributionOfTaxa.txt";
open IN, "<$taxFile" or die("Can't open $taxFile for reading.\n");
my $tax = {};

while ( my $line = <IN> ) {
   if ( $line =~ /^\S/ ) {
      chomp $line;
      (my $taxType, my $taxVal, my $prefixAndSpecies, my $dirName) = split(/\t/, $line);
      if ( $taxType eq 'phylum' or $taxType eq 'class' or $taxType eq 'order' or $taxType eq 'family' or $taxType eq 'genus' ) {
	 $tax->{$dirName}->{$taxType} = $taxVal;
      }
   }
}

close IN;

$date = `date`;
chomp $date;
print "done.  [$date]\n\n";

###########################################################################################
# 
# Query the phylogenetic info of input reads by BLASTN
#
###########################################################################################
print "Querying taxa in BLAST for input reads [$date]...\n\n";
print LOG "Querying taxa in BLAST for input reads [$date]...\n\n";

my $blOut="$Bin/rawBlastOut_$input_file";
system("blastn -query $input -out $blOut -db $blastDB -num_threads 2 -outfmt 7");

my $inBlastScore = {};
my $bestHitScore = {};
my $bestHitName ={};
my $allUniqTaxon={};

open BLASTIN, "<$blOut" or die("Can't open $blOut for reading.\n");

while ( my $line = <BLASTIN> ) {
	if ( $line !~ /^#/ ) {
		chomp $line;
		my @fields = split(/\t/, $line);
		my $queryID = $fields[0];
		my $matchName = $fields[1];
		$matchName =~ s/\.\d+$//;
      		my $currentScore = $fields[10];
		if ( not $inBlastScore->{$queryID}->{$matchName} or $inBlastScore->{$queryID}->{$matchName} > $currentScore ) {
	 		$inBlastScore->{$queryID}->{$matchName} = $currentScore;
   		}
		if ( not $bestHitScore->{$queryID} or $bestHitScore->{$queryID} > $currentScore ) {
			$bestHitScore->{$queryID} = $currentScore;
			$bestHitName->{$queryID} = $matchName;
		}
   	} 	
}

close BLASTIN;

#my $inResultFile="$Bin/result_blast_all_$input_file";
#open OUT, ">$inResultFile" or die("Can't open $inResultFile for writing.\n");
#print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";
#foreach my $queryID ( sort { $a cmp $b } keys %$inBlastScore ) {
#	foreach my $spName (keys %{$inBlastScore->{$queryID}}) {
#		print OUT "$queryID\t$spName\t$inBlastScore->{$queryID}->{$spName}\t$tax->{$spName}->{'genus'}\t$tax->{$spName}->{'family'}\t$tax->{$spName}->{'order'}\t$tax->{$spName}->{'class'}\t$tax->{$spName}->{'phylum'}\n";
#	}
#}
#close OUT;

my $bestHitResultFile="$Bin/result_blast_besthit_$input_file";
open OUT, ">$bestHitResultFile" or die("Can't open $bestHitResultFile for writing.\n");
print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";

foreach my $queryID (keys %$bestHitScore) {
	my $spName=$bestHitName->{$queryID};
	print OUT "$queryID\t$spName\t$bestHitScore->{$queryID}\t$tax->{$spName}->{'genus'}\t$tax->{$spName}->{'family'}\t$tax->{$spName}->{'order'}\t$tax->{$spName}->{'class'}\t$tax->{$spName}->{'phylum'}\n";
}

close OUT;

push(@tmpFiles,$blOut);

# Find out the dominant bins based on BLAST results
my $numRdAllBin=0;
my $allSpecCnt={};
my $numRdTaxonAllBin={};
foreach my $id (keys %$bestHitName) {
	my $name=$bestHitName->{$id};

	if(not $numRdTaxonAllBin->{'species'}->{$name}) {
		$numRdTaxonAllBin->{'species'}->{$name}=1;
	} else {
		$numRdTaxonAllBin->{'species'}->{$name}++;
	}

	my $genus=$tax->{$name}->{'genus'};
	if(not $numRdTaxonAllBin->{'genus'}->{$genus}) {
		$numRdTaxonAllBin->{'genus'}->{$genus}=1;
	} else {
		$numRdTaxonAllBin->{'genus'}->{$genus}++;
	}

	my $family=$tax->{$name}->{'family'};
	if(not $numRdTaxonAllBin->{'family'}->{$family}) {
		$numRdTaxonAllBin->{'family'}->{$family}=1;
	} else {
		$numRdTaxonAllBin->{'family'}->{$family}++;
	}

	my $order=$tax->{$name}->{'order'};
	if(not $numRdTaxonAllBin->{'order'}->{$order}) {
		$numRdTaxonAllBin->{'order'}->{$order}=1;
	} else {
		$numRdTaxonAllBin->{'order'}->{$order}++;
	}

	my $class=$tax->{$name}->{'class'};
	if(not $numRdTaxonAllBin->{'class'}->{$class}) {
		$numRdTaxonAllBin->{'class'}->{$class}=1;
	} else {
		$numRdTaxonAllBin->{'class'}->{$class}++;
	}

	my $phylum=$tax->{$name}->{'phylum'};
	if(not $numRdTaxonAllBin->{'phylum'}->{$phylum}) {
		$numRdTaxonAllBin->{'phylum'}->{$phylum}=1;
	} else {
		$numRdTaxonAllBin->{'phylum'}->{$phylum}++;
	}

	$numRdAllBin++;
	$allSpecCnt->{$name}=$numRdTaxonAllBin->{'species'}->{$name};
}

$date = `date`;
chomp $date;
print "done.  [$date]\n";
print LOG "done.  [$date]\n";

###########################################################################################
#
# Call BLASTN to further classify the bins based on e-value
#
###########################################################################################
$date = `date`;
chomp $date;
print "Querying taxa in BLAST for GPD outputs [$date]...\n\n";
print LOG "Querying taxa in BLAST for GPD outputs [$date]...\n\n";

# Dominant bins, in the type of file_name=>bin_name
my $domBins={};
# all the scores and matched name in the whole genome
my $resultSummary="$Bin/result_summary_$output";

open SUMMARY, ">$resultSummary" or die("Can't open $resultSummary for writing.\n");

print "$realNumofBins bins have been identified.\n\n";
print SUMMARY "$realNumofBins bins have been identified.\n\n";

my $binCnt=0;
foreach my $blinput (@gpdOutputs) {

	my $blastFile="$Bin/.rawBlast_$output.$binCnt";
	system("blastn -query $blinput -out $blastFile -db $blastDB -num_threads 2 -outfmt 7");

	# Read the raw BLAST E-values.
	my $blastScore = {};
	my $blastMatch = {};

	open BLASTIN, "<$blastFile" or die("Can't open $blastFile for reading.\n");

	while ( my $line = <BLASTIN> ) {
		if ( $line !~ /^#/ ) {
			chomp $line;
			my @fields = split(/\t/, $line);
			my $queryID = $fields[0];
			my $matchName = $fields[1];
			$matchName =~ s/\.\d+$//;
      			my $currentScore = $fields[10];
			if ( not $blastScore->{$queryID} or $blastScore->{$queryID} > $currentScore ) {
	 			$blastScore->{$queryID} = $currentScore;
				$blastMatch->{$queryID} = $matchName;
      			}
   		} 	
	}

	close BLASTIN;

	my $resultFile="$Bin/result_blast_$output.$binCnt";
	open OUT, ">$resultFile" or die("Can't open $resultFile for writing.\n");
	print OUT "QUERY_ID\tBEST_MATCH\tSCORE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";
	foreach my $queryID ( sort { $a cmp $b } keys %$blastScore ) {
		my $spName=$blastMatch->{$queryID};
	   	print OUT "$queryID\t$spName\t$blastScore->{$queryID}\t$tax->{$spName}->{'genus'}\t$tax->{$spName}->{'family'}\t$tax->{$spName}->{'order'}\t$tax->{$spName}->{'class'}\t$tax->{$spName}->{'phylum'}\n";
	}
	close OUT;

	# Calculate p-value
	#
	# 			annotated with term A	Annotated w/o term A
	# Reference gene set		m			n-m
	# Target gene set		M			N-M
	#
	# n: total # of reads in this bin
	# m: # of reads of species A
	# M: # of reads of species A in all bins
	# N: total # of reads in all bins
	my $numRdThisBin=0;
	my $taxonCnt={};
	my $speciesCnt={};
	foreach my $id (keys %$blastMatch) {
		my $name=$blastMatch->{$id};

		if(not $taxonCnt->{'species'}->{$name}) {
			$taxonCnt->{'species'}->{$name}=1;
		} else {
			$taxonCnt->{'species'}->{$name}++;
		}

		my $genus=$tax->{$name}->{'genus'};
		if(not $taxonCnt->{'genus'}->{$genus}) {
			$taxonCnt->{'genus'}->{$genus}=1;
		} else {
			$taxonCnt->{'genus'}->{$genus}++;
		}

		my $family=$tax->{$name}->{'family'};
		if(not $taxonCnt->{'family'}->{$family}) {
			$taxonCnt->{'family'}->{$family}=1;
		} else {
			$taxonCnt->{'family'}->{$family}++;
		}

		my $order=$tax->{$name}->{'order'};
		if(not $taxonCnt->{'order'}->{$order}) {
			$taxonCnt->{'order'}->{$order}=1;
		} else {
			$taxonCnt->{'order'}->{$order}++;
		}

		my $class=$tax->{$name}->{'class'};
		if(not $taxonCnt->{'class'}->{$class}) {
			$taxonCnt->{'class'}->{$class}=1;
		} else {
			$taxonCnt->{'class'}->{$class}++;
		}

		my $phylum=$tax->{$name}->{'phylum'};
		if(not $taxonCnt->{'phylum'}->{$phylum}) {
			$taxonCnt->{'phylum'}->{$phylum}=1;
		} else {
			$taxonCnt->{'phylum'}->{$phylum}++;
		}

		$speciesCnt->{$name}=$taxonCnt->{'species'}->{$name};
		$numRdThisBin++;
	}

	die "ERROR: failed to find executable phyper under $Bin!\n" if !-x "$Bin/phyper";
	print SUMMARY "\n===============================\n";

	my $pValue = {};
	foreach my $id (keys %$blastMatch) {
		my $name=$blastMatch->{$id};
		my $genus=$tax->{$name}->{'genus'};
		my $family=$tax->{$name}->{'family'};
		my $order=$tax->{$name}->{'order'};
		my $class=$tax->{$name}->{'class'};
		my $phylum=$tax->{$name}->{'phylum'};

		my $numOtherSpecs=$numRdAllBin-$numRdTaxonAllBin->{'species'}->{$name};
		$pValue->{$name}->{'species'}=`$Bin/phyper -k $taxonCnt->{'species'}->{$name} -n1 $numRdTaxonAllBin->{'species'}->{$name} -n2 $numOtherSpecs -t $numRdThisBin`;
		chomp $pValue->{$name}->{'species'};
		$pValue->{$name}->{'species'}=1-$pValue->{$name}->{'species'};

		my $numOtherGenus=$numRdAllBin-$numRdTaxonAllBin->{'genus'}->{$genus};
		$pValue->{$name}->{'genus'}=`$Bin/phyper -k $taxonCnt->{'genus'}->{$genus} -n1 $numRdTaxonAllBin->{'genus'}->{$genus} -n2 $numOtherGenus -t $numRdThisBin`;
		chomp $pValue->{$name}->{'genus'};
		$pValue->{$name}->{'genus'}=1-$pValue->{$name}->{'genus'};

		my $numOtherFamily=$numRdAllBin-$numRdTaxonAllBin->{'family'}->{$family};
		$pValue->{$name}->{'family'}=`$Bin/phyper -k $taxonCnt->{'family'}->{$family} -n1 $numRdTaxonAllBin->{'family'}->{$family} -n2 $numOtherFamily -t $numRdThisBin`;
		chomp $pValue->{$name}->{'family'};
		$pValue->{$name}->{'family'}=1-$pValue->{$name}->{'family'};

		my $numOtherOrder=$numRdAllBin-$numRdTaxonAllBin->{'order'}->{$order};
		$pValue->{$name}->{'order'}=`$Bin/phyper -k $taxonCnt->{'order'}->{$order} -n1 $numRdTaxonAllBin->{'order'}->{$order} -n2 $numOtherOrder -t $numRdThisBin`;
		chomp $pValue->{$name}->{'order'};
		$pValue->{$name}->{'order'}=1-$pValue->{$name}->{'order'};

		my $numOtherClass=$numRdAllBin-$numRdTaxonAllBin->{'class'}->{$class};
		$pValue->{$name}->{'class'}=`$Bin/phyper -k $taxonCnt->{'class'}->{$class} -n1 $numRdTaxonAllBin->{'class'}->{$class} -n2 $numOtherClass -t $numRdThisBin`;
		chomp $pValue->{$name}->{'class'};
		$pValue->{$name}->{'class'}=1-$pValue->{$name}->{'class'};

		my $numOtherPhylum=$numRdAllBin-$numRdTaxonAllBin->{'phylum'}->{$phylum};
		$pValue->{$name}->{'phylum'}=`$Bin/phyper -k $taxonCnt->{'phylum'}->{$phylum} -n1 $numRdTaxonAllBin->{'phylum'}->{$phylum} -n2 $numOtherPhylum -t $numRdThisBin`;
		chomp $pValue->{$name}->{'phylum'};
		$pValue->{$name}->{'phylum'}=1-$pValue->{$name}->{'phylum'};
	}

	# sort speciesCnt by value in descending order
	print "GPD BIN $binCnt ($numRdThisBin reads, totally $numRdAllBin reads in all bins):\n\n";
	print SUMMARY "GPD BIN $binCnt ($numRdThisBin reads, totally $numRdAllBin reads in all bins):\n\n";
	print SUMMARY "SPECIES\tABUNDANCE(p-value)\tGENUS(p-value)\tFAMILY(p-value)\tORDER(p-value)\tCLASS(p-value)\tPHYLUM(p-value)\n";
	foreach my $name ( sort {$speciesCnt->{$b}<=>$speciesCnt->{$a}} keys %$speciesCnt) { 
		my $genus=$tax->{$name}->{'genus'};
		my $family=$tax->{$name}->{'family'};
		my $order=$tax->{$name}->{'order'};
		my $class=$tax->{$name}->{'class'};
		my $phylum=$tax->{$name}->{'phylum'};

		#	print "$name\t$taxonCnt->{'species'}->{$name}($pValue->{$name}->{'species'})\t$taxonCnt->{'genus'}->{$genus}($pValue->{$name}->{'genus'})\t$taxonCnt->{'family'}->{$family}($pValue->{$name}->{'family'})\t$taxonCnt->{'order'}->{$order}($pValue->{$name}->{'order'})\t$taxonCnt->{'class'}->{$class}($pValue->{$name}->{'class'})\t$taxonCnt->{'phylum'}->{$phylum}($pValue->{$name}->{'phylum'})\n"; 
		print SUMMARY "$name\t$taxonCnt->{'species'}->{$name}($pValue->{$name}->{'species'})\t$taxonCnt->{'genus'}->{$genus}($pValue->{$name}->{'genus'})\t$taxonCnt->{'family'}->{$family}($pValue->{$name}->{'family'})\t$taxonCnt->{'order'}->{$order}($pValue->{$name}->{'order'})\t$taxonCnt->{'class'}->{$class}($pValue->{$name}->{'class'})\t$taxonCnt->{'phylum'}->{$phylum}($pValue->{$name}->{'phylum'})\n"; 
	}	
	
	push(@tmpFiles,$blastFile); # record temp files

	$binCnt++;
}

print SUMMARY "THE WHOLE GENOME:\n\n";
print SUMMARY "BEST_MATCH\tABUNDANCE\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";
foreach my $name ( sort {$allSpecCnt->{$b}<=>$allSpecCnt->{$a}} keys %$allSpecCnt) { 
	my $genus=$tax->{$name}->{'genus'};
	my $family=$tax->{$name}->{'family'};
	my $order=$tax->{$name}->{'order'};
	my $class=$tax->{$name}->{'class'};
	my $phylum=$tax->{$name}->{'phylum'};

	print SUMMARY "$name\t$allSpecCnt->{$name}\t$genus($numRdTaxonAllBin->{'genus'}->{$genus})\t$family($numRdTaxonAllBin->{'family'}->{$family})\t$order($numRdTaxonAllBin->{'order'}->{$order})\t$class($numRdTaxonAllBin->{'class'}->{$class})\t$phylum($numRdTaxonAllBin->{'phylum'}->{$phylum})\n";
}	

close SUMMARY;

$date = `date`;
chomp $date;
print "done.  [$date]\n";
print LOG "done.  [$date]\n";

close LOG;

# delete temp files
foreach my $file (@tmpFiles) {
	system("rm -f $file");
}

#####################################################
#	SUBROUTINES
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
