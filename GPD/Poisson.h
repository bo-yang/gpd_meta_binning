/*
 * Poisson.h
 *
 *  Created on: Oct 5, 2012
 *      Author: bo
 */

#ifndef POISSON_H_
#define POISSON_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <ctype.h>
#include <zlib.h>
#include <cmath>
#include <limits>
#include <time.h>
#include "kseq.h"

using namespace std;

/*
 * Mixture of Poisson
 */
class Poisson {
public:
	int* bin;	// result of binning for each x_i

	// functions
	Poisson(int poisnSize);
	virtual ~Poisson();
	void EMinit();
	int EMrun();
	int LoadData(char* inputFile);	// Read gene sequence data and do statistics
	int LoadPoissonData(char* inputFile);	// Read Poisson data
	void Classify();
	void WriteResult(char* resultFile);
	void InitAlpha(); // Initialize alpha
	void InitLambda();	// initialize lambda

private:
	static const int DEFAULT_POISSON_SIZE;
	static const double CONVERGE_THRESHHOLD;
	static const int KMER_SIZE;	// a fixed size of k-tuple
	static const int KMER_SIZE_MAX;
	static const double VERY_SMALL_NUMBER;
	int* x;	// occurrences of k-mers
	long kmer_size;	 // number of different k-mers, n
	int poisn_size;	// number of Poisson distributions, m
	double* alpha;	// the probability that a sample is drawn from one of the m Poisson distributions
	double* lambda;	// the parameters for the m Poisson distributions
	double* initLambda; // initial lambda
	double** z;	// P(y_j=j|X,Theta)
	map<string,int> kmer_occur;
	int lambda_upperlimit;
	int* subtotal;	// subtotal of x_i for each Poisson distribution

	double CalcProb(int i, int j);	// calculate probability
	double Calc_z(int i, int j);	// calculate z_ij
	void ShowLambda();
	void Merge(int* input, long p, long r); // merge sort
	void Merge_sort(int* input, long p, long r); // merge sort
};

#endif /* POISSON_H_ */
