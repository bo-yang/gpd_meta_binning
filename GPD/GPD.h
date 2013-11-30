/*
 * GPD.h
 *
 * Generalized Poisson Distribution.
 *
 *  Created on: Aug 23, 2012
 *      Author: boyang
 */

#ifndef GENPOISS_H_
#define GENPOISS_H_

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
#include <assert.h>
#include <vector>
#include <gsl/gsl_cdf.h>
#include "kseq.h"

using namespace std;

struct KmerBin
{
	int bin;	// the bin #
	vector<int> kmers;	// kmer occurs belongs to this bin
	vector<double> probs;	// probability of kmer occurrence
	int bin_size;	// # of different kmers in this bin
	int kmer_size;	// # total kmers in this bin
};

class GPD
{
public:
	// variables
	static const int PREDEFINED_GPD_SIZE;
	static const double CONVERGE_THRESHHOLD; // threshhold of EM convergence
	static const int KMER_SIZE_MAX;	// maximum number of simulated k-mer occurrences
	static const int LAMBDA_SIZE;	// number of lambdas, default 2 for each GPD.
	static const int KMER_LEN;	// a fixed size of k-mer
	static const double NEWTON_EPS;	// convergence criteria for Newton's method
	static const int EM_MAX_ITER; // max number of EM iteration.
	static const int NEWTON_MAX_ITER; // max number of Newton's method iteration.
	static const int READ_COVERAGE;	// read coverage of a sequence
	static const double VERY_SMALL_NUMBER;
	static const int FILE_NAME_LEN;
	static const double VALID_KMER_PERCENT;	// threshold of valid k-mer percentage
	bool recordKmerLen;

	// functions
	GPD(int gpdSize=PREDEFINED_GPD_SIZE,bool simdata=false,int kmerlen=KMER_LEN,double em_eps=CONVERGE_THRESHHOLD,double kmer_percent=VALID_KMER_PERCENT);
	~GPD();
	void EMinit();
	int EMrun();
	int LoadData(char* inputFile);	// Read gene sequence data and do statistics
	int LoadPoissonData(char* inputFile);	// Read Poisson data
	void WriteResult_kmer(char* resultFile);
	void WriteResult_read(char* inputFile, char* resultFile);
	void WriteScore(char* resultFile);
	void InitAlpha();	// initialize alpha
	void InitAlpha(double* vector);	// initialize alpha
	void InitLambda();	// initialize lambda
	void InitLambda(double** matrix);	// initialize lambda
	void Binning(char* inputFile);	// assign reads to bins
	void SetGPDsize(int num);	// set the number of different GPDs
	void SetKmerPercent(double percent);	// set the valid kmer percentage
	void ReleaseMemory();	// Delete all memory allocated dynamically
	void Classify_kmer();	// assign k_mer to bins
	double CalcCHI();	// calculate Calinski-Harabasz Index
	double CalcSWC();	// calculate Silhouette Width Criterion
	double Likelihood(); // Calculate the likelihood of L(X,Y,L|Theta)
	double LogLikelihood(); // Calculate the log likelihood of L(X,Y,L|Theta)
	double BIC(double lh); // Calculate the Bayesian Information Criterion
	double BIC(double L1, double L2, int m1, int m2, int n);	// calculate the Bayesian Information criterion
	double IndexI();	// Index I score
	double XBI();	// Xie Beni Index
	static double phyper(unsigned int k, unsigned int n1, unsigned int n2, unsigned int t);	// cummulative hypergeometric function

private:
	// variables
	int gpd_size;	// number of different GPDs(m).
	int uniq_kmer_size;	// number of different k-mers(n).
	int all_kmer_size;	// number of all k-mers(N).
	vector<long> x;		// observed data, x_i is the unique occurrence of the i-th k-tuple in all reads.
	vector<long> l;	// observed data, occurrences of each x_i
	map<string,int> kmer_occur;	// k-mer and occurrences
	double* alpha;	// the probability that a sample x_i is from a GPD.
	double* initAlpha; // initial value of alpha
	double** lambda;	// parameters of GPD.
	double** initLambda;	// initial value of lambda, it useful to record the workable value of lambda.
	double** z;	// expectation of hidden value y_ij
	map<int,int> binOfKmer;	// result of binning for each kmer x_i. First: x[i]; second: bin of x[i]
	vector<KmerBin> kmerbin;	// results of each bin
	bool isSimData;	// flag of simulated data
	int kmer_len;	// the length of k continuous sequence letters
	double em_thresh; // Converge threshhold of EM algorithm
	double valid_kmer_percent; // percentage of valid kmers, to determine the value of initial lambda

	// functions
	double Calc_z(unsigned int i,unsigned int j);	// calculate z
	double CalcProb(unsigned int i,unsigned int k);	// calculate probability
	double prob(int x,unsigned int j);	// calculate probability by specifying x
	double cdf(unsigned int x,unsigned int j);	// cumulative distribution
	double Calc_w(int j);	// parameter of f(lambda_j2)
	double f_lambda_j2(int j);
	double df_lambda_j2(int j);
	int NewtonMethod(int lambda_index); // Newton's method to calculate lambda_j2
	int CheckParam(int j); // check if the updated lambdas meet the requirement: 0<lambda_j2<1
	void ShowLambda();
	void Merge(int* input, long p, long r); // merge sort
	void Merge_sort(int* input, long p, long r); // merge sort
	int GetMaxBin(int* a, int len);
	int GetMaxBin(double* a, int len);
	double dist(int x1,double p1,int x2,double p2);	// distance between points (x1,p1) and (x2,p2)
	double dist(double x,double y);	// distance between vectors x and y
	double CalcJSD(double* x,double* y,int len);	// calculate Jensen-Shannon Divergence(JSD)
	double CalcKLD(double* x,double* y,int len);	// calculate Kullback-Leibler Divergence(KLD)
	double CalcJSD(double x,double y);	// calculate Jensen-Shannon Divergence(JSD)
	double CalcKLD(double x,double y);	// calculate Kullback-Leibler Divergence(KLD)
	void TuneGPD();	// insert 0 & 1 into x[i] so that k-mer occurs fit GPD
	double CalcReadProb(vector<int> &kmers, int j);	// calculate the probability that a read belongs to a bin
	int GetIndexofx(int occur);	// Get the index of kmer in x
};

#endif /* GENPOISS_H_ */
