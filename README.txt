This is the source code of metagenomic binning. Metagenomics is the study of
microbial communities sampled directly from their natural environment, without
prior culturing. The goal this project is to cluster metagenomic reads into as many species as possible. 

The major assumption is that the mixed metagenomic reads are in different
generalized Poisson distributions(GPD). Therefore, a GPD model is built and
the parameters are estimated by EM algorithm(codes
under GPD dir). Details of the model can be found in the documents.

Explantion of major codes/tools:
1. gpd_reads
Tool to generate synthetic data in Generalized Poisson Distribution:
	./tools/gpd_reads -i ./Genome/NC_010693.fna -o gpdReads_3_genomes_random_cover1to6to12_small.fna -s 1 -c 1.2

Command to build gpd_reads(GSL is required):
	g++ -Wall -o gpd_reads gpd_reads.cpp -lgsl -lgslcblas -lm

2. ./tools/phyper
Tool to calculate the hypergeometric distribution of Poisson distribution. 

3. ./GPD/gpd
	g++ -g -o gpd *.h *.cpp -lz -lgsl -lgslcblas -lm

