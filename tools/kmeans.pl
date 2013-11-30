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
$output="$basename" if (!$output);

my @tempFile;

# Generate k-mer vectors for each read
print "Counting $kmerlen-mer occurrences....\n";
my $matfile="$basename\_$kmerlen-mer.mat";
system("./gen_kmer_mat.pl -i $input -o $matfile -k $kmerlen");
push(@tempFile,$matfile);
my $kmers=ReadMat($matfile); # Vector of kmer occurrences

# Call GraphLab k-means method and find the optimum clustering
print "Clustering reads by K-means....\n";
my $clst=2;
my $cnt=0;
my $bestscore=-99999;
my $optclust=$clst;
my $optclustfile;
while($cnt<3) {
	print "\t$clst clusters => ";
	my $clustfile="$output\_$clst\_clusts_cluster";
	my $datafile="$output\_$clst\_clusts_data";
	system("./kmeans --data=$matfile --clusters=$clst --id=1 --output-clusters=$clustfile --output-data=$datafile >/dev/null 2>&1");
	$datafile=`ls $datafile*`; # it maybe renamed as MetaCluster_A2_out.0_2_clusts_data_1_of_1
	chomp($datafile);
	push(@tempFile,$clustfile);
	push(@tempFile,$datafile);

	my $cents=ReadCents($clustfile); # Read centroids
	my $clusts=ReadData($datafile); # Clustered read IDs

	my $scoreII=II($clst,$cents,$kmers,$clusts);
	my $scoreXB=XB($clst,$cents,$kmers,$clusts);
	print "Index I: $scoreII; Xie Beni: $scoreXB\n";
	if($bestscore<=$scoreII) {
		$bestscore=$scoreII;
		$optclust=$clst;
		$optclustfile=$datafile;
		$cnt=0;
	} else {
		# If the new scores are smaller than the best score for continuous 3 times,
		# stop the iteration.
		$cnt++;
	}
	
	$clst++;
}

# Calculate p-value
print "Parsing clustering results and calculating p-value....\n";
system("./pvalue_kmeans_withid.pl -t taxonomyData -i $input -o $output.$optclust-clusters.pvalue -c $optclustfile 1>/dev/null");

# Remove temporary files
# foreach my $file (@temFile) {
# 	system("rm -f $file");
# }

print "Done! $optclust clusters found.\n";

#########################################################
#	SUBROUTINES
#########################################################

#
# Index I score
#
sub II {
	my ($K,$cents,$kmers,$clusts)=@_;
	my $p=2;

	my $E1=0;
	foreach my $id (@{$clusts->{1}}) {
		if(length($kmers->{$id})==0){
			$id="0$id"; # TEMPORARY, the ID 07 would be 7, which would lead to an error.
			#print "II: $id\n"; # TEST ONLY
		}
		$E1 += DistL2($kmers->{$id},$cents->{1});
	}

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

	#print "K=$K, E1=$E1, EK=$EK, DK=$DK\n"; # TEST ONLY
	my $score=(1/$K*$E1/$EK*$DK)**$p;

	return $score;
}

#
# Xie Beni score
#
sub XB {
	my ($K,$cents,$kmers,$clusts)=@_;
	my $N=keys %$kmers; # total number of points
	my $numer=0;
	foreach my $clst (sort keys %$clusts) {
		foreach my $id (@{$clusts->{$clst}}) {
			if(length($kmers->{$id})==0){
				$id="0$id"; # TEMPORARY, the ID 07 would be 7, which would lead to an error.
				#print "II: $id\n"; # TEST ONLY
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

#
# L1 norm
#
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

#
# L2 norm
#
sub DistL2 {
	my ($p1_ref,$p2_ref)=@_;
	my @p1=split(/ /,$p1_ref);
	my @p2=split(/ /,$p2_ref);

	if(scalar(@p1)!=scalar(@p2)) {
		print "ERROR: DistL2: dimensions of two points not match!\n";
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
# 	Input: cluster file containing centroids
# 	Output: hash(cluster) of hashes(centroid) 
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
# 	Output: hash(cluster) of arrays(read IDs)
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
