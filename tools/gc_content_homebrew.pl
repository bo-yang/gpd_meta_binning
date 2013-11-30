#!/usr/bin/perl

###
# GC_content.pl - This tool counts GC content.
#
###

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my $input;
my $output;

GetOptions('input|i=s'  => \$input,
	   'output|o=s' => \$output
   );

die("ERROR: no input reads specied!\n") if (!$input);
my $basename=`basename $input`;
chomp($basename);
$output="$basename.gc" if (!$output);

open IN, "<$input" or die("ERROR: cannot open file $input!");
open OUT, ">$output" or die("ERROR: cannot create file $output!\n");

my %nt=('A',0, # Neucleotides
	'T',0,
	'C',0,
	'G',0);

# Simply count the GC content
#while (my $line = <IN>) {
#	chomp $line;
#	if($line =~ /^>/) {
#		# Record read ID
#		(my $read,my $loc,my $gi, my $gi_num, my $ref, my $ref_id, my $source_1)=split(/\|/,$line);
#		$source_1 =~ s/ /_/g;
#		$source_1 =~ s/,.*//g; # remove comma
#		$source_1 =~ s/^_//g;  # leading _
#		$source_1 =~ s/_chromosome.*//g;
#		$source_1 =~ s/_plasmid.*//g;
#		$source_1 =~ s/_incision.*//g;
#		my $spec=$source_1;
#
#		print OUT "$read\t$spec\t";
#
#	} else {
#		%nt=( 	# Neucleotides
#			'A',0, 
#			'T',0,
#			'C',0,
#			'G',0
#		);
#
#		for(my $i=0;$i<length($line)-1;++$i) {
#			my $char=substr($line,$i,1);
#			$nt{$char}++;		
#		}
#
#		my $gc_cont=($nt{'G'}+$nt{'C'})/($nt{'A'}+$nt{'T'}+$nt{'G'}+$nt{'C'})*100;
#		print OUT "$gc_cont\n";
#	} # end of if
#}

# Get the distribution of GC content per line
my $num=0;
my $gc_dist={};
while (my $line = <IN>) {
	chomp $line;
	if($line =~ /^>/) {
		# Record read ID
		(my $read_id,my $loc,my $gi, my $gi_num, my $ref, my $ref_id, my $source_1)=split(/\|/,$line);
		$source_1 =~ s/ /_/g;
		$source_1 =~ s/,.*//g; # remove comma
		$source_1 =~ s/^_//g;  # leading _
		$source_1 =~ s/_chromosome.*//g;
		$source_1 =~ s/_plasmid.*//g;
		$source_1 =~ s/_incision.*//g;
		my $spec=$source_1;

	} else {
		%nt=( 	# Neucleotides
			'A',0, 
			'T',0,
			'C',0,
			'G',0
		);

		for(my $i=0;$i<length($line)-1;++$i) {
			my $char=substr($line,$i,1);
			$nt{$char}++;		
		}

		my $gc_cont=($nt{'G'}+$nt{'C'})/($nt{'A'}+$nt{'T'}+$nt{'G'}+$nt{'C'})*100;
		if( not $gc_dist->{$gc_cont}) {
			$gc_dist->{$gc_cont}=1;
		} else {
			$gc_dist->{$gc_cont}++;
		}
		$num++;
	} # end of if
}

foreach my $gc_cont (sort {$gc_dist->{$b} <=> $gc_dist->{$a}} keys %{$gc_dist}) {
	my $ratio=$gc_dist->{$gc_cont}/$num;
	print OUT "$gc_cont\t$ratio\n";
}

close IN;
close OUT;
