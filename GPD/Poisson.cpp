/*
 * Poisson.cpp
 *
 *  Created on: Oct 5, 2012
 *      Author: bo
 */

#include "Poisson.h"

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

const int Poisson::DEFAULT_POISSON_SIZE=3;
const double Poisson::CONVERGE_THRESHHOLD=0.0000001;
const int Poisson::KMER_SIZE=12;
const int Poisson::KMER_SIZE_MAX=50000;
const double Poisson::VERY_SMALL_NUMBER=numeric_limits<double>::min();

Poisson::Poisson(int poisnSize=DEFAULT_POISSON_SIZE) {
	x=NULL;
	z=NULL;
	alpha=NULL;
	lambda=NULL;
	initLambda=NULL;
	bin=NULL;
	kmer_size=0;
	poisn_size=poisnSize;
	lambda_upperlimit=20;
	subtotal=NULL;
}

Poisson::~Poisson() {
	if(x!=NULL)
		delete [] x;

	if(z!=NULL)
	{
		for(int i=0;i<kmer_size;++i)
			delete [] z[i];
		delete [] z;
	}

	if(alpha!=NULL)
		delete [] alpha;

	if(lambda!=NULL)
		delete [] lambda;
	if(initLambda!=NULL)
		delete [] initLambda;
	if(subtotal!=NULL)
		delete [] subtotal;

	x=NULL;
	z=NULL;
	alpha=NULL;
	lambda=NULL;
	subtotal=NULL;
}

// Read FASTA data and do some statistics
int Poisson::LoadData(char* inputFile)
{
	gzFile fp;
	kseq_t *seq; // genome sequence
	int l=0;

	fp = gzopen(inputFile, "r"); // STEP 2: open the file handler
	seq = kseq_init(fp); // STEP 3: initialize seq

	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence

		// count the occurrences of k-tuples in each pseudo genome
		for(unsigned int i=0;i<strlen(seq->seq.s)-KMER_SIZE;++i)
		{
			string k_mer(seq->seq.s+i,KMER_SIZE);
			pair<map<string,int>::iterator,bool> ret;
			ret=kmer_occur.insert(pair<string,int>(k_mer,1));
			if (ret.second==false)
			{
				// k_mer already exists in kmer_occur
				++kmer_occur[k_mer];
			}
		} // end of for
	}
	kseq_destroy(seq); // STEP 5: destroy seq
	gzclose(fp); // STEP 6: close the file handler

	// Allocate memory for x
	kmer_size=kmer_occur.size();
	x=new int[kmer_size];
	memset(x,0,sizeof(int)*kmer_size);

	long i=0;
	map<string,int>::iterator it;
	for ( it=kmer_occur.begin() ; it != kmer_occur.end(); it++ )
	{
		if(lambda_upperlimit<(*it).second)
			lambda_upperlimit=(*it).second;
		if((*it).second>1)		// exclude k-mers occurring once
			x[i++]=(*it).second;
	}
	kmer_size=i;  // update the size of k-mers after exluding 1s.

	/*	fstream fs;
	fs.open("poisson_data", fstream::out);
	for(int i=0;i<kmer_size;++i)
	{
		fs<<x[i]<<" ";
		if(i!=0 && i%20==0)
			fs<<endl;
	}
	fs<<endl;
	fs.close();
	 */
	Merge_sort(x,0,kmer_size-1);

	return kmer_size;
}

// read simulated poisson data
int Poisson::LoadPoissonData(char* inputFile)
{
	kmer_size=0;

	// Allocate memory for x
	x=new int[KMER_SIZE_MAX];
	memset(x,0,sizeof(int)*KMER_SIZE_MAX);

	fstream fs;
	fs.open(inputFile, fstream::in);
	if(fs.is_open() == false)
	{
		cout<<"ERROR: failed to open file "<<inputFile<<"!\n";
		return -1;
	}

	// Read sequences
	while(!fs.eof())
	{
		fs>>x[kmer_size];
		if(lambda_upperlimit<x[kmer_size])
			lambda_upperlimit=x[kmer_size];
		kmer_size++;
	} // end of while

	kmer_size--;	// remove the count of EOF

	fs.close();

	Merge_sort(x,0,kmer_size-1);

	return kmer_size;
}

void Poisson::Merge(int* input, long p, long r)
{
    long mid = floor((p + r) / 2);
    long i1 = 0;
    long i2 = p;
    long i3 = mid + 1;

    // Temp array
    int* temp=new int[r-p+1];

    // Merge in sorted form the 2 arrays
    while ( i2 <= mid && i3 <= r )
        if ( input[i2] < input[i3] )
            temp[i1++] = input[i2++];
        else
            temp[i1++] = input[i3++];

    // Merge the remaining elements in left array
    while ( i2 <= mid )
        temp[i1++] = input[i2++];

    // Merge the remaining elements in right array
    while ( i3 <= r )
        temp[i1++] = input[i3++];

    // Move from temp array to master array
    for ( int i = p; i <= r; i++ )
        input[i] = temp[i-p];

    delete [] temp;
}

// inputs:
//	p - the start index of array input
//	r - the end index of array input
void Poisson::Merge_sort(int* input, long p, long r)
{
    if ( p < r )
    {
        long mid = floor((p + r) / 2);
        Merge_sort(input, p, mid);
        Merge_sort(input, mid + 1, r);
        Merge(input, p, r);
    }
}

void Poisson::EMinit()
{
	// initialize EM algorithm
	// X, kmer_size, gpd_size should be already configured in LoadData()

	// initialization of alpha
	alpha=new double[poisn_size];

	// allocate memory for lambda,y and z
	lambda=new double[poisn_size];
	initLambda=new double[poisn_size];

	z=new double*[kmer_size];
	for(int i=0;i<kmer_size;++i)
	{
		z[i]=new double[poisn_size];
		memset(z[i],0,sizeof(double)*poisn_size);
	}

	// allocate memory
	bin=new int[kmer_size];
	subtotal=new int[poisn_size];
}

void Poisson::InitAlpha()
{
//	srand(time(NULL));
	for(int i=0;i<poisn_size;++i)
	{
		alpha[i]=(double)1/poisn_size;
//		alpha[i]=(double)rand()/(double)RAND_MAX/(1+VERY_SMALL_NUMBER)/(double)poisn_size+VERY_SMALL_NUMBER;
	}
}

// initialize lambda.
// Randomly choose several k-mer occurrences, calculate the mean of
// randomly-chosen data, and assign it to lambda.
void Poisson::InitLambda()
{
	InitAlpha();

	/* initialize random seed: */
	srand(time(NULL));

	int i=0;
	for(int j=0;j<poisn_size;++j)
	{
		int numSample=kmer_size*alpha[j];
		double x_mean=0.0;
		for(int cnt=0;cnt<numSample;++cnt)
		{
			x_mean+=x[i+cnt];
		}
		i+=numSample;
		lambda[j]=x_mean/(double)numSample;
		initLambda[j]=lambda[j];
	}
}
/*
int Poisson::EMrun()
{
	int em_ret=1;

	double* oldLambda=new double[poisn_size];

	// TODO: test which is better: update parameters separately or in a bunch?
	// update alpha & lambda
	while(true)
	{
		double numer=0.0;
		double denom=0.0;
		// update z
		for(int i=0;i<kmer_size;++i)
			for(int k=0;k<poisn_size;++k)
				z[i][k]=Calc_z(i,k);

		for(int j=0;j<poisn_size;++j)
		{
			// record the value of lambda
			oldLambda[j]=lambda[j];

			double subtotal_z=0.0;
			double numerator=0.0;
			for(int i=0;i<kmer_size;++i)
			{
				numerator += z[i][j]*x[i];
				subtotal_z += z[i][j];
			}
			alpha[j]=subtotal_z/(double)kmer_size;
			lambda[j]=numerator/subtotal_z;

			numer += pow(lambda[j]-oldLambda[j],2);
			denom += pow(lambda[j],2);
		} // end of for
		// check for convergence
		//			double oldLambda_j=oldLambda[j];
		//			double lambda_j=lambda[j];
		double rate=sqrt(numer)/sqrt(denom);
		if(rate<CONVERGE_THRESHHOLD)
		{
			//				ShowLambda();	// TEST ONLY
			break;
		}
	} // end of while

	delete [] oldLambda;

	return em_ret;
}
*/

int Poisson::EMrun()
{
	int em_ret=1;

	double* oldLambda=new double[poisn_size];

	// TODO: test which is better: update parameters separately or in a bunch?
	// update alpha & lambda
	for(int j=0;j<poisn_size;++j)
	{
		while(true)
		{
			// record the value of lambda
			oldLambda[j]=lambda[j];

			// update z
			for(int i=0;i<kmer_size;++i)
				for(int k=0;k<poisn_size;++k)
					z[i][k]=Calc_z(i,k);

			double subtotal_z=0.0;
			double numerator=0.0;
			for(int i=0;i<kmer_size;++i)
			{
				numerator += z[i][j]*x[i];
				subtotal_z += z[i][j];
			}
			alpha[j]=subtotal_z/(double)kmer_size;
			lambda[j]=numerator/subtotal_z;

			// check for convergence
			//			double oldLambda_j=oldLambda[j];
			//			double lambda_j=lambda[j];
			double rate=abs((lambda[j]-oldLambda[j])/lambda[j]);
			if(rate<CONVERGE_THRESHHOLD)
			{
				//				ShowLambda();	// TEST ONLY
				break;
			}
		} // end of while
	} // end of for

	delete [] oldLambda;

	return em_ret;
}

// inputs:
//		i - index of x
//		j - index of lambda and alpha
double Poisson::CalcProb(int i, int j)
{
	double log_xi=0.0;
	for(int k=1;k<=x[i];++k)
		log_xi += log(k);

	double prob=exp(x[i]*log(lambda[j])-lambda[j]-log_xi);
	if(isnan(prob)||isinf(prob))
	{
//		cout<<"ERROR: data out of range in Poisson::CalcProb()."<<endl;
		prob=VERY_SMALL_NUMBER;
	}

	return prob;
}

double Poisson::Calc_z(int i, int j)
{
	double denom=0.0;
	for(int k=0;k<poisn_size;++k)
		denom += alpha[k]*CalcProb(i,k);

	double frac=alpha[j]*CalcProb(i,j)/denom;
	if(isnan(frac)||isinf(frac))
	{
		cout<<"ERROR: data out of range in Poisson::Calc_z()."<<endl;
	}

	return frac;
}

void Poisson::ShowLambda()
{
	for(int j=0;j<poisn_size;++j)
		cout<<lambda[j]<<" ";
	cout<<endl;
}

void Poisson::Classify()
{
	memset(subtotal,0,sizeof(int)*poisn_size);
	// Get the latest P(y=j)
	for(int i=0;i<kmer_size;++i)
	{
		z[i][0]=Calc_z(i,0);
		bin[i]=0;
		double max=z[i][0];
		for(int j=1;j<poisn_size;++j)
		{
			z[i][j]=Calc_z(i,j);
			if(max<z[i][j])
			{
				max=z[i][j];
				bin[i]=j;
			}
		}
		subtotal[bin[i]]++;
	} // end of for
}

void Poisson::WriteResult(char* resultFile)
{
	fstream fs;
	fs.open(resultFile, fstream::out|fstream::app);
	if(fs.is_open() == false)
	{
		cout<<"ERROR: failed to open file "<<resultFile<<"!\n";
		return;
	}

	fs<<"\n================================================\n";
	int final_bins=poisn_size;
	for(int j=0;j<poisn_size;++j)
	{
		if(subtotal[j]==0)
			final_bins--;
	}
	fs<<"Classified species: "<<final_bins<<endl;
	fs<<"Initial Lambdas:"<<endl;
	for(int j=0;j<poisn_size;++j)
		fs<<initLambda[j]<<" ";
	fs<<endl<<endl;
	fs<<"Final Lambdas:"<<endl;
	for(int j=0;j<poisn_size;++j)
		fs<<lambda[j]<<" ";
	fs<<endl<<endl;

	fs<<"Final Alphas:"<<endl;
	for(int j=0;j<poisn_size;++j)
		fs<<alpha[j]<<" ";
	fs<<endl<<endl;

	for(int j=0;j<poisn_size;++j)
	{
		fs<<"Poisson distribution "<<j<<"(total "<<subtotal[j]<<"):"<<endl;
/*		for(int i=0;i<kmer_size;++i)
		{
			if(bin[i]==j)
				fs<<x[i]<<" ";
		}*/
		fs<<endl;
	}

	fs.close();
}
