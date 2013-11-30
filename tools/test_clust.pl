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
use MCE;

my $numclust=3;
my $datafile;
my $clustfile;
my $matfile;
my $output;

GetOptions('data|d=s' => \$datafile,
	   'cent|c=s' => \$clustfile,
	   'mat|m=s'  => \$matfile,
	   'k=i'      => \$numclust,
	   'output|o=s'	      => \$output
   );

my $cents=ReadCents($clustfile);
my $kmers=ReadMat($matfile);
my $clusts=ReadData($datafile);

my $scoreII=II($numclust,$cents,$kmers,$clusts);
my $scoreXB=XB($numclust,$cents,$kmers,$clusts);
print "Index I: $scoreII\n";
print "Xie Beni: $scoreXB\n";

#########################################################
#	SUBROUTINES
#########################################################

# Index I score
sub II {
	my ($K,$cents,$kmers,$clusts)=@_;
	my $p=2;

	my $E1=0;
	my $pm = Parallel::ForkManager->new(5);
	# data structure retrieval and handling
	$pm -> run_on_finish ( # called BEFORE the first call to start()
		sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
	  
			# retrieve data structure from child
			if (defined($data_ref)) {  # children are not forced to send anything
				$E1+=$data_ref;  # child passed a string reference
	      		} else {  # problems occuring during storage or retrieval will throw a warning
	        		print "ERROR: got no data from subprocess!\n";
	      		}
	    	}
	);

	foreach my $id (@{$clusts->{1}}) {
		$pm->start() and next;

		if(length($kmers->{$id})==0){
			$id="0$id"; # TEMPORARY
			#print "II: $id\n"; # TEST ONLY
		}
		#$E1 += DistL2($kmers->{$id},$cents->{1});
		my $sum+=DistL2($kmers->{$id},$cents->{1});
		# send it back to the parent process
    		$pm->finish(0, \$sum);  # note that it's a scalar REFERENCE, not the scalar itself
	}
	$pm->wait_all_children;

	my $EK=0;
	foreach my $clst (sort keys %$clusts) {
		foreach my $id (@{$clusts->{$clst}}) {
			if(length($kmers->{$id})==0){
				$id="0$id"; # TEMPORARY
				#print "II: $id\n"; # TEST ONLY
			}

			$EK += DistL2($kmers->{$id},$cents->{$clst});
		}
	}

	my $DK=0;
	foreach my $i (sort keys %$cents) {
		foreach my $j (sort keys %$cents) {
			my $dist=DistL2($cents->{$i},$cents->{$j});
			if($DK<$dist){
				$DK=$dist;
			}
		}
	}

	print "K=$K, E1=$E1, EK=$EK, DK=$DK\n"; # TEST ONLY
	my $score=(1/$K*$E1/$EK*$DK)**$p;

	return $score;
}

# Xie Beni score
sub XB {
	my ($K,$cents,$kmers,$clusts)=@_;
	my $N=keys %$kmers; # total number of points
	my $numer=0;
	foreach my $clst (sort keys %$clusts) {
		foreach my $id (@{$clusts->{$clst}}) {
			if(length($kmers->{$id})==0){
				$id="0$id"; # TEMPORARY
				#print "XB: $id\n"; # TEST ONLY
			}

			$numer+=DistL2($kmers->{$id},$cents->{$clst})**2;
		}
	}

	my $denom=999999999999;
	foreach my $i (sort keys %$cents) {
		foreach my $j (sort keys %$cents) {
			my $dist=DistL2($cents->{$i},$cents->{$j})**2;
			if($i!=$j && $denom>$dist){
				$denom=$dist;
			}
		}
	}

	my $score=$numer/($N*$denom);
	return $score;
}

# L1 norm
sub DistL1 {
	my ($p1_ref,$p2_ref)=@_;
	my @p1=split(/ /,$p1_ref);
	my @p2=split(/ /,$p2_ref);

	if(scalar(@p1)!=scalar(@p2)) {
		print "ERROR: DistL1: dimensions of two points not match!\n";
		return 'inf';
	}

	my $dist=0;
	for(my $i=0;$i<scalar(@p1);$i++) {
		$dist+=abs($p1[$i]-$p2[$i]);
	}

	return $dist;
}

# L2 norm
sub DistL2 {
	my ($p1_ref,$p2_ref)=@_;
	my @p1=split(/ /,$p1_ref);
	my @p2=split(/ /,$p2_ref);

	if(scalar(@p1)!=scalar(@p2)) {
		print "ERROR: DistL2: dimensions of two points not match!\n";
		#print "Point 1:\n@p1\n";
		#print "Point 2:\n@p2\n";
		return 'inf';
	}

	my $dist=0;
	for(my $i=0;$i<scalar(@p1);$i++) {
		$dist += ($p1[$i]-$p2[$i])**2;
	}
	$dist=sqrt($dist);

	return $dist;
}

# Read centroids from clust file
sub ReadCents {
	my $clustfile=$_[0];
	my $cents={};

	open IN, "<$clustfile" or die("Can't open $clustfile for reading.\n");
	while ( my $line = <IN> ) {
		chomp $line;
		(my $clst,my $centroid)=split(/\t/,$line);
		$cents->{$clst}=$centroid;
	}
	close IN;
	
	return($cents);
}

# Load read vectors of k-mer occurrences
# 	Input: file of read vectors
# 	Output: hash(read IDs) of hashes(k-mers)
sub ReadMat {
	my $matfile=$_[0];
	my $kmers={};

	open IN, "<$matfile" or die("Can't open $matfile for reading.\n");
	while ( my $line = <IN> ) {
		chomp $line;
		(my $id,my $vector)=split(/\t/,$line);
		$kmers->{$id}=$vector;
	}
	close IN;

	return($kmers);
}

# Load clustered data file
# 	Input: data file
# 	Output: hash(cluster) of hashes(read IDs)
sub ReadData {
	my $datafile=$_[0];
	my $clusts={};

	open IN, "<$datafile" or die("Can't open $datafile for reading.\n");
	while ( my $line = <IN> ) {
		chomp $line;
		(my $id,my $clst)=split(/\t/,$line);
		push(@{$clusts->{$clst}},$id);
	}
	close IN;

	return($clusts);
}
