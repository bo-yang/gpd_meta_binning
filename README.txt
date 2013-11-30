1. gpd_reads
Tool to generate synthetic data in Generalized Poisson Distribution:
	./tools/gpd_reads -i ./Genome/NC_010693.fna -o gpdReads_3_genomes_random_cover1to6to12_small.fna -s 1 -c 1.2

Command to build gpd_reads(GSL is required):
	g++ -Wall -o gpd_reads gpd_reads.cpp -lgsl -lgslcblas -lm

2. ./tools/phyper
Tool to calculate the hypergeometric distribution of Poisson distribution. 

3. gpd
	g++ -g -o gpd *.h *.cpp -lz -lgsl -lgslcblas -lm

4. gpdbl.pl

