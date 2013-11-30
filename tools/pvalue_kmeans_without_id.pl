#!/usr/bin/perl

###
# pvalue_kmeans.pl - This tool analyse the reads clustered by K-means.
#
# Usage: pvalue_kmeans.pl -i <input_reads> -c <clustered_reads> [-o <out_file>] 
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';
use Math::GSL::CDF qw/:all/;

my $input;
my $output;
my $clustfile;
my $taxPath;
GetOptions('input|d=s'	 => \$input,
	   'output|o=s'  => \$output,
   	   'clustfile|c=s'=> \$clustfile,
	   'taxonomy|t=s' => \$taxPath);

die("ERROR: no input input specied!\n") if (!$input);

open IN, "<$input" or die("ERROR: cannot open file $input!");

# Handle signals
my $CleanupDone   = 0;
$SIG{'ABRT'} = 'doCleanup';

# Ignore following signals.
$SIG{'HUP'}  = 'doCleanup';
$SIG{'INT'}  = 'doCleanup';
$SIG{'QUIT'} = 'doCleanup';

# p-value threshold
my $EPS=0.05;

###########################################################################################
# 
# Read organism names so we can generate results files comparable to what BLAST will output.
# 
###########################################################################################

print "Loading taxonomy...\n";

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

###########################################################################################
# 
# Reading cluster file.
# 
###########################################################################################

print "Reading cluster file...\n";

open IN, "<$clustfile" or die("Can't open $clustfile for reading.\n");
my $clusts={};
my $cluststats={};
my $cnt=1;
while ( my $line = <IN> ) {
	chomp $line;
	(my $clst, my $junk)=split(/\./,$line);
	$clst =~ s/^\s+//; # trim leading spaces
	$clusts->{$cnt}=$clst;
	if ( not $cluststats->{$clst} ) {
		$cluststats->{$clst}=1;
	} else {
		$cluststats->{$clst}++;
	}
	$cnt++;
}

close IN;

my $numclust=keys %$cluststats; # number of different clusters

###########################################################################################
# 
# Count species and get the p-values
#
###########################################################################################
print "Calculating p-values...\n";

my $clustcladestats={}; # Store clades info of each cluster
my $cladestats={}; # Store all clades info
my $clustspecstats={}; # Store the species info
my $specstats={}; # Store unique species info

open IN, "<$input" or die("ERROR: cannot open file $input!");

$cnt=1;
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

		if (not $clustcladestats->{$clusts->{$cnt}}->{$spec}) {
			$clustcladestats->{$clusts->{$cnt}}->{$spec}=1;
		} else {
			$clustcladestats->{$clusts->{$cnt}}->{$spec}++;
		}
		if (not $cladestats->{$spec}) {
			$cladestats->{$spec}=1;
		} else {
			$cladestats->{$spec}++;
		}


		my $genus=$tax->{$spec}->{'genus'};
		if (not $clustcladestats->{$clusts->{$cnt}}->{$genus}) {
			$clustcladestats->{$clusts->{$cnt}}->{$genus}=1;
		} else {
			$clustcladestats->{$clusts->{$cnt}}->{$genus}++;
		}
		if (not $cladestats->{$genus}) {
			$cladestats->{$genus}=1;
		} else {
			$cladestats->{$genus}++;
		}

		my $family=$tax->{$spec}->{'family'};
		if (not $clustcladestats->{$clusts->{$cnt}}->{$family}) {
			$clustcladestats->{$clusts->{$cnt}}->{$family}=1;
		} else {
			$clustcladestats->{$clusts->{$cnt}}->{$family}++;
		}
		if (not $cladestats->{$family}) {
			$cladestats->{$family}=1;
		} else {
			$cladestats->{$family}++;
		}

		my $order=$tax->{$spec}->{'order'};
		if (not $clustcladestats->{$clusts->{$cnt}}->{$order}) {
			$clustcladestats->{$clusts->{$cnt}}->{$order}=1;
		} else {
			$clustcladestats->{$clusts->{$cnt}}->{$order}++;
		}
		if (not $cladestats->{$order}) {
			$cladestats->{$order}=1;
		} else {
			$cladestats->{$order}++;
		}

		my $class=$tax->{$spec}->{'class'};
		if (not $clustcladestats->{$clusts->{$cnt}}->{$class}) {
			$clustcladestats->{$clusts->{$cnt}}->{$class}=1;
		} else {
			$clustcladestats->{$clusts->{$cnt}}->{$class}++;
		}
		if (not $cladestats->{$class}) {
			$cladestats->{$class}=1;
		} else {
			$cladestats->{$class}++;
		}

		my $phylum=$tax->{$spec}->{'phylum'};
		if (not $clustcladestats->{$clusts->{$cnt}}->{$phylum}) {
			$clustcladestats->{$clusts->{$cnt}}->{$phylum}=1;
		} else {
			$clustcladestats->{$clusts->{$cnt}}->{$phylum}++;
		}
		if (not $cladestats->{$phylum}) {
			$cladestats->{$phylum}=1;
		} else {
			$cladestats->{$phylum}++;
		}

		$clustspecstats->{$clusts->{$cnt}}->{$spec}=$clustcladestats->{$clusts->{$cnt}}->{$spec};
		$specstats->{$spec}=$cladestats->{$spec};

		$cnt++;
	}

}
close IN;
my $totreads=$cnt;

my $basename=`basename $input`;
chomp($basename);
$output="$basename.$numclust\_clusters.pvalue" if (!$output);

open OUT, ">$output" or die("ERROR: cannot create file $output!\n");

my $pvalue={};
my $cladeerr={}; # wrong classified reads for all clades
my $flagclst={};
my $err={};
# Calculate p-value
foreach my $clst (sort keys %$cluststats) {
	print OUT "===============================\n";
	print OUT "Cluster $clst ($cluststats->{$clst} reads, totally $totreads reads in all clusters):\n\n";
	print OUT "SPECIES\tABUNDANCE(p-value)\tGENUS(p-value)\tFAMILY(p-value)\tORDER(p-value)\tCLASS(p-value)\tPHYLUM(p-value)\n";
	my $thisclust=$clustcladestats->{$clst};
	foreach my $spec (sort {$clustspecstats->{$clst}->{$b}<=>$clustspecstats->{$clst}->{$a}} keys %{$clustspecstats->{$clst}}) {
		my $genus=$tax->{$spec}->{'genus'};
		my $family=$tax->{$spec}->{'family'};
		my $order=$tax->{$spec}->{'order'};
		my $class=$tax->{$spec}->{'class'};
		my $phylum=$tax->{$spec}->{'phylum'};

		$pvalue->{$clst}->{$spec}=1-gsl_cdf_hypergeometric_P($thisclust->{$spec},$cladestats->{$spec},$totreads-$cladestats->{$spec},$cluststats->{$clst});
		if($pvalue->{$clst}->{$spec} > $EPS && $flagclst->{$spec}!=$clst) {
			$cladeerr->{$spec}+=$thisclust->{$spec};
			$flagclst->{$spec}=$clst; # species are unique, so only add once.
		}

		$pvalue->{$clst}->{$genus}=1-gsl_cdf_hypergeometric_P($thisclust->{$genus},$cladestats->{$genus},$totreads-$cladestats->{$genus},$cluststats->{$clst});
		if($pvalue->{$clst}->{$genus} > $EPS && $flagclst->{$genus}!=$clst) {
			$cladeerr->{$genus}+=$thisclust->{$genus};
			$flagclst->{$genus}=$clst; # only add once
		}

		$pvalue->{$clst}->{$family}=1-gsl_cdf_hypergeometric_P($thisclust->{$family},$cladestats->{$family},$totreads-$cladestats->{$family},$cluststats->{$clst});
		if($pvalue->{$clst}->{$family} > $EPS && $flagclst->{$family}!=$clst) {
			$cladeerr->{$family}+=$thisclust->{$family};
			$flagclst->{$family}=$clst;  # only add once
		}

		$pvalue->{$clst}->{$order}=1-gsl_cdf_hypergeometric_P($thisclust->{$order},$cladestats->{$order},$totreads-$cladestats->{$order},$cluststats->{$clst});
		if($pvalue->{$clst}->{$order} > $EPS && $flagclst->{$order}!=$clst) {
			$cladeerr->{$order}+=$thisclust->{$order};
			$flagclst->{$order}=$clst;  # only add once
		}

		$pvalue->{$clst}->{$class}=1-gsl_cdf_hypergeometric_P($thisclust->{$class},$cladestats->{$class},$totreads-$cladestats->{$class},$cluststats->{$clst});
		if($pvalue->{$clst}->{$class} > $EPS && $flagclst->{$class}!=$clst) {
			$cladeerr->{$class}+=$thisclust->{$class};
			$flagclst->{$class}=$clst;  # only add once
		}

		$pvalue->{$clst}->{$phylum}=1-gsl_cdf_hypergeometric_P($thisclust->{$phylum},$cladestats->{$phylum},$totreads-$cladestats->{$phylum},$cluststats->{$clst});
		if($pvalue->{$clst}->{$phylum} > $EPS && $flagclst->{$phylum}!=$clst) {
			$cladeerr->{$phylum}+=$thisclust->{$phylum};
			$flagclst->{$phylum}=$clst; # only add once
		}

		print OUT "$spec\t$thisclust->{$spec}($pvalue->{$clst}->{$spec})\t$thisclust->{$genus}($pvalue->{$clst}->{$genus})\t$thisclust->{$family}($pvalue->{$clst}->{$family})\t$thisclust->{$order}($pvalue->{$clst}->{$order})\t$thisclust->{$class}($pvalue->{$clst}->{$class})\t$thisclust->{$phylum}($pvalue->{$clst}->{$phylum})\n";
	}

}

###########################################################################################
# 
# Calculate clustering errors.
#
###########################################################################################
print "Calculating errors...\n";

print OUT "\n===============================\n";
print OUT "THE WHOLE GENOME:\n\n";
print OUT "SPECIES(err)\tABUNDANCE\tGENUS(err)\tFAMILY(err)\tORDER(err)\tCLASS(err)\tPHYLUM(err)\n";
foreach my $spec ( sort {$specstats->{$b}<=>$specstats->{$a}} keys %{$specstats}) { 
	my $genus=$tax->{$spec}->{'genus'};
	my $family=$tax->{$spec}->{'family'};
	my $order=$tax->{$spec}->{'order'};
	my $class=$tax->{$spec}->{'class'};
	my $phylum=$tax->{$spec}->{'phylum'};

	$err->{'spec'} += $cladeerr->{$spec};
	$err->{'genus'} += $cladeerr->{$genus};
	$err->{'family'} += $cladeerr->{$family};
	$err->{'order'} += $cladeerr->{$order};
	$err->{'class'} += $cladeerr->{$class};
	$err->{'phylum'} += $cladeerr->{$phylum};

	$cladeerr->{$spec}/=$cladestats->{$spec};
	$cladeerr->{$genus}/=$cladestats->{$genus};
	$cladeerr->{$family}/=$cladestats->{$family};
	$cladeerr->{$order}/=$cladestats->{$order};
	$cladeerr->{$class}/=$cladestats->{$class};
	$cladeerr->{$phylum}/=$cladestats->{$phylum};

	my $spec_err=sprintf("%.4f",$cladeerr->{$spec});
	my $genus_err=sprintf("%.4f",$cladeerr->{$genus});
	my $family_err=sprintf("%.4f",$cladeerr->{$family});
	my $order_err=sprintf("%.4f",$cladeerr->{$order});
	my $class_err=sprintf("%.4f",$cladeerr->{$class});
	my $phylum_err=sprintf("%.4f",$cladeerr->{$phylum});

	print OUT "$spec($spec_err)\t$specstats->{$spec}\t$genus($genus_err)\t$family($family_err)\t$order($order_err)\t$class($class_err)\t$phylum($phylum_err)\n"; 
}

# errors in general
$err->{'spec'} /= $totreads;
$err->{'genus'} /= $totreads;
$err->{'family'} /= $totreads;
$err->{'order'} /= $totreads;
$err->{'class'} /= $totreads;
$err->{'phylum'} /= $totreads;

my $spec_err=sprintf("%.4f",$err->{'spec'});
my $genus_err=sprintf("%.4f",$err->{'genus'});
my $family_err=sprintf("%.4f",$err->{'family'});
my $order_err=sprintf("%.4f",$err->{'order'});
my $class_err=sprintf("%.4f",$err->{'class'});
my $phylum_err=sprintf("%.4f",$err->{'phylum'});
print OUT "\n===============================\n";
print OUT "GENERAL ERRORS:\n\n";
print OUT "SPECIES\tGENUS\tFAMILY\tORDER\tCLASS\tPHYLUM\n";
print OUT "$spec_err\t$genus_err\t$family_err\t$order_err\t$class_err\t$phylum_err\n";

close OUT;

print "Done!\n";

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
