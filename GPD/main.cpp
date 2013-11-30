/*
 * main.cpp
 *
 *  Created on: Sep 9, 2012
 *      Author: boyang
 */
#include <iostream>
#include "GPD.h"
#include "Poisson.h"

using namespace std;

const int MAX_EM_ITER=1000;
const double KMER_PERCENT_INCR=0.05;

// calculate the Bayesian Information criterion:
//	BIC=2*ln(L1/L2)-(m1-m2)*ln(n)
double BIC(double L1, double L2, int m1, int m2, int n)
{
	double bic=2*log(L1/L2)-(m1-m2)*log(n);
	return bic;
}

int main(int argc, char* argv[])
{
	char* inputFile=NULL;
	char* outputFile=NULL;
	int composition_len=GPD::KMER_LEN;
	int bin_size=GPD::PREDEFINED_GPD_SIZE;
	bool loadPoissonData=false;
	bool isGPD=true;
	bool autoBinning=false;
	bool recordKmerLen=false;
	double em_eps=GPD::CONVERGE_THRESHHOLD;
	double kmer_percent=GPD::VALID_KMER_PERCENT;

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-input") == 0)
		{
			i++;
			inputFile = argv[i];
		}
		else if (strcmp(argv[i], "-kmer_len") == 0)
		{
			i++;
			composition_len = atoi(argv[i]);
			recordKmerLen=true;
		}
		else if (strcmp(argv[i], "-output") == 0)
		{
			i++;
			outputFile = argv[i];
		}
		else if (strcmp(argv[i], "-bin_size") == 0)
		{
			i++;
			bin_size = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-lp") == 0)
		{
			loadPoissonData=true;
		}
		else if (strcmp(argv[i], "-poisson") == 0)
		{
			isGPD=false;
		}
		else if (strcmp(argv[i], "-em_threshhold") == 0)
		{
			i++;
			em_eps=atof(argv[i]);
		}
		else if (strcmp(argv[i], "-kmer_percent") == 0)
		{
			i++;
			kmer_percent=atof(argv[i]);
		}
		else if (strcmp(argv[i], "-auto") == 0)
		{
			autoBinning=true;
		}
	}

	if(inputFile==NULL)
	{
		cout<<"ERROR: no input file was specified!"<<endl;
		exit(1);
	}
	if(inputFile==NULL)
	{
		cout<<"ERROR: no output file was specified!"<<endl;
		exit(1);
	}

	if(autoBinning)
	{
		bin_size=2; // starting number of bins

		GPD gp(bin_size,loadPoissonData,composition_len,em_eps,kmer_percent);
		gp.recordKmerLen=recordKmerLen;
		cout<<"Start binning now, please wait...."<<endl;
		if(loadPoissonData)
			gp.LoadPoissonData(inputFile);
		else
			gp.LoadData(inputFile);

		int numBelowOpt=0;	// number of run times below optimum score
		double bestScore=-1;
		int optBinSize=0;	// optimum number of bins
		while(numBelowOpt<3){
			// prepare for EM algorithm
			int ret=0;
			int iter=0;
			kmer_percent=GPD::VALID_KMER_PERCENT;
			while(!ret && iter<MAX_EM_ITER)
			{
				gp.ReleaseMemory();	// release memory
				gp.SetKmerPercent(kmer_percent);
				// initialize parameters
				gp.EMinit();
				// run EM algorithm
				ret=gp.EMrun();
				iter++;
				kmer_percent+=KMER_PERCENT_INCR;
			}

			if(iter == MAX_EM_ITER)
			{
				gp.ReleaseMemory();	// release memory
				break;
			}

			gp.Classify_kmer();
			gp.WriteScore(outputFile); //// TEST ONLY ////
			double scoreCHI=gp.CalcCHI();	// Temporarily only use CHI
			//double scoreII=gp.IndexI(); // Temporarily only use Index I
			double score=scoreCHI; //// TEST ONLY ////
			if(bestScore<score)
			{
				bestScore=score;
				optBinSize=bin_size;
			} else {
				numBelowOpt++;
			}

			gp.ReleaseMemory();	// release memory
			gp.SetGPDsize(++bin_size);
		}

		// binning based on optimized parameters
		gp.ReleaseMemory();	// release memory
		gp.SetGPDsize(optBinSize);
		kmer_percent=GPD::VALID_KMER_PERCENT;
		gp.SetKmerPercent(kmer_percent);

		int iter=0;
		int ret=0;
		while(!ret)
		{
			// initialize parameters
			gp.EMinit();
			// run EM algorithm
			ret=gp.EMrun();

			if(iter>2*MAX_EM_ITER)
			{
				cout<<"FATAL ERROR: "<<inputFile<<" cannot be classified with k-mer="<<composition_len<<"."<<endl;
				break;
			}
			iter++;
		}

		// Write binning results
		if(loadPoissonData)
		{
			gp.WriteResult_kmer(outputFile);
		} else {
			gp.WriteResult_read(inputFile,outputFile);
		}
		cout<<optBinSize<<" species have been identified."<<endl;

	} else {
		if(isGPD)
		{
			GPD gp(bin_size,loadPoissonData,composition_len,em_eps,kmer_percent);
			gp.recordKmerLen=recordKmerLen;
			cout<<"Loading data now, please wait...."<<endl;
			if(loadPoissonData)
				gp.LoadPoissonData(inputFile);
			else
				gp.LoadData(inputFile);

			cout<<"Running EM-Algorithm now, please wait...."<<endl;
			// prepare for EM algorithm
			int ret=0;
			int iter=0;
			while(!ret)
			{
				cout<<"EM-Algorithm: loop "<<iter<<endl;
				// initialize parameters
				gp.EMinit();
				// run EM algorithm
				ret=gp.EMrun();
				if(iter>2*MAX_EM_ITER)
				{
					cout<<"FATAL ERROR: "<<inputFile<<" cannot be classified with k-mer="<<composition_len<<"."<<endl;
					return 1;
				}
				iter++;
			}
			// binning based on optimized parameters
			if(loadPoissonData)
			{
				gp.WriteResult_kmer(outputFile);
			} else {
				gp.WriteResult_read(inputFile,outputFile);
			}
			cout<<"Done!"<<endl;
		} else {
			Poisson gp(bin_size);
			cout<<"Loading data now, please wait...."<<endl;
			if(loadPoissonData)
				gp.LoadPoissonData(inputFile);
			else
				gp.LoadData(inputFile);

			cout<<"Running EM-Algorithm now, please wait...."<<endl;
			// prepare for EM algorithm
			gp.EMinit();
			int ret=0;
			while(!ret)
			{
				//			cout<<"EM-Algorithm: loop "<<em_cnt++<<endl;
				// initialize parameters
				gp.InitLambda();
				// run EM algorithm
				ret=gp.EMrun();
			}
			// binning based on optimized parameters
			gp.Classify();
			gp.WriteResult(outputFile);
			cout<<"Done!"<<endl;
		}
	}

	return 0;
}
