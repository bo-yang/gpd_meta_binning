/*
 * GenPoiss.cpp
 *
 *  Created on: Aug 23, 2012
 *      Author: boyang
 */

#include "GPD.h"

const int GPD::PREDEFINED_GPD_SIZE=3;
const int GPD::LAMBDA_SIZE=2;
const int GPD::KMER_LEN=16;
const int GPD::KMER_SIZE_MAX=50000; // naively set the maximum number of samples
const double GPD::CONVERGE_THRESHHOLD=0.08;
const double GPD::NEWTON_EPS=0.0001;
const int GPD::EM_MAX_ITER=5000;
const int GPD::NEWTON_MAX_ITER=5000;
const int GPD::READ_COVERAGE=3;
const double GPD::VERY_SMALL_NUMBER=numeric_limits<double>::min();
const int GPD::FILE_NAME_LEN=256;
const double GPD::VALID_KMER_PERCENT=0.01;	// 0.12 for 3 species

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

GPD::GPD(int gpdSize,bool simdata,int kmerlen,double em_eps,double kmer_percent)
{
	// initialization
	gpd_size=gpdSize;
	uniq_kmer_size=0;
	all_kmer_size=0;
	z=NULL;
	lambda=NULL;
	alpha=NULL;
	initAlpha=NULL;
	initLambda=NULL;
	isSimData=simdata;
	kmer_len=kmerlen;
	recordKmerLen=false;
	em_thresh=em_eps;
	valid_kmer_percent=kmer_percent;

	srand(time(NULL));
}

GPD::~GPD()
{
	ReleaseMemory();
}

void GPD::ReleaseMemory()
{
	// release dynamically allocated memory of lambda
	if(lambda!=NULL)
	{
		for(int i=0;i<gpd_size;++i)
		{
			delete [] lambda[i];
			delete [] initLambda[i];
		}
		delete [] lambda;
		delete [] initLambda;
	}

	// release dynamically allocated memory of alpha
	if(alpha!=NULL)
	{
		delete [] alpha;
	}
	if(initAlpha!=NULL)
	{
		delete [] initAlpha;
	}

	if(z!=NULL)
	{
		for(int i=0;i<uniq_kmer_size;++i)
			delete [] z[i];
		delete [] z;
	}

	lambda=NULL;
	alpha=NULL;
	initAlpha=NULL;
	z=NULL;
	initLambda=NULL;

	binOfKmer.clear();
	kmerbin.clear();
}

void GPD::SetGPDsize(int num)
{
	gpd_size=num;
}

void GPD::SetKmerPercent(double percent)
{
	valid_kmer_percent=percent;
}

// Read FASTA data and do some statistics
int GPD::LoadData(char* inputFile)
{
	gzFile fp;
	kseq_t *seq; // genome sequence
	int len=0;

	fp = gzopen(inputFile, "r"); // STEP 2: open the file handler
	seq = kseq_init(fp); // STEP 3: initialize seq

	while ((len = kseq_read(seq)) >= 0) { // STEP 4: read sequence
		// count the occurrences of k-tuples in each pseudo genome
		for(unsigned int i=0;i<strlen(seq->seq.s)-kmer_len;++i)
		{
			string k_mer(seq->seq.s+i,kmer_len);
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

	map<int,int> x_stat;	// first number: unique number in x; second number: total number of each unique x
	map<string,int>::iterator it;
	for ( it=kmer_occur.begin() ; it != kmer_occur.end(); it++ )
	{
		if((*it).second>1)		// exclude k-mers occurring once
		{
			pair<map<int,int>::iterator,bool> ret;
			ret=x_stat.insert(pair<int,int>((*it).second,1));
			if (ret.second==false)
			{
				// k_mer already exists in kmer_occur
				++x_stat[(*it).second];
			}
		}
	}

	uniq_kmer_size=x_stat.size();
	all_kmer_size=0;

	int cnt=0;	// reserve positions for 0 & 1
	map<int,int>::iterator itx;
	for(itx=x_stat.begin();itx!=x_stat.end();++itx)
	{
		// elements in x should be sorted ascendingly when inserting elements into x_stat.
		x.insert(x.end(),(*itx).first);
		l.insert(l.end(),(*itx).second);
		all_kmer_size+=l[cnt];
		cnt++;
	}
	//kmer_bin_size=x[cnt-1]+1;	// index binOfKmer with k-mer occurs, i.e. binOfKmer[x[i]].

	fstream fs;
	fs.open("kmer_occur_stat", fstream::out);
	map<int,int>::iterator its;
	fs<<all_kmer_size<<" kmers have been loaded."<<endl;
	fs<<"kmer_occur \t stats"<<endl;
	for (int i=0;i<uniq_kmer_size;++i)
	{
		fs<<x[i]<<": \t "<<l[i]<<endl;
	}
	fs<<endl;
	fs.close();

	return uniq_kmer_size;
}

// Load simulated data
int GPD::LoadPoissonData(char* inputFile)
{
	uniq_kmer_size=0;

	fstream fs;
	fs.open(inputFile, fstream::in);
	if(fs.is_open() == false)
	{
		cout<<"ERROR: failed to open file "<<inputFile<<"!\n";
		return -1;
	}

//	int tmp;
	map<int,int> x_stat;	// first number: unique number in x; second number: total number of each unique x
	// Read sequences
	/*	while(!fs.eof())
	{
		fs>>tmp;
		pair<map<int,int>::iterator,bool> ret;
		ret=x_stat.insert(pair<int,int>(tmp,1));
		if (ret.second==false)
		{
			// k_mer already exists in kmer_occur
			++x_stat[tmp];
		}
	} // end of while
	 */
	int cnt=0;
	while(!fs.eof())
	{
		int tmp;
		fs>>tmp;
		x.insert(x.end(),tmp);
		fs>>tmp;
		l.insert(l.end(),tmp);
		all_kmer_size+=l[cnt];
		cnt++;
	}
	uniq_kmer_size=cnt-1;
	fs.close();

	//	uniq_kmer_size=x_stat.size();
	//	all_kmer_size=0;

	/*	int cnt=0;
	map<int,int>::iterator itx;
	for(itx=x_stat.begin();itx!=x_stat.end();++itx)
	{
		// elements in x should be sorted ascendingly when inserting elements into x_stat.
		if((*itx).second>0)
		{
			x.insert(x.end(),(*itx).first);
			l.insert(l.end(),(*itx).second);
			all_kmer_size+=l[cnt];
			cnt++;
		}
	}
	 */
	fs.open("kmer_occur_stat", fstream::out);
	map<int,int>::iterator its;
	fs<<"kmer_occur \t stats"<<endl;
	for (int i=0;i<uniq_kmer_size;++i)
	{
		fs<<x[i]<<": \t "<<l[i]<<endl;
		if(x[i]<=0){
			cout<<"ERROR: invalid input data!"<<endl;
		}
	}
	fs<<endl;
	fs.close();

	return uniq_kmer_size;
}

// Preparation for EM algorithm
void GPD::EMinit()
{
	// initialize EM algorithm
	// X, uniq_kmer_size, gpd_size should be already configured in LoadData()

	// initialization of alpha
	if(alpha==NULL)
		alpha=new double[gpd_size];

	// allocate memory for lambda,y and z
	if(lambda==NULL)
	{
		lambda=new double*[gpd_size];
		for(int i=0;i<gpd_size;++i)
		{
			lambda[i]=new double[LAMBDA_SIZE];
		}
	}
	if(initLambda==NULL)
	{
		initLambda=new double*[gpd_size];
		for(int i=0;i<gpd_size;++i)
		{
			initLambda[i]=new double[LAMBDA_SIZE];
		}
	}

	if(z==NULL)
	{
		z=new double*[uniq_kmer_size+2];
		for(int i=0;i<uniq_kmer_size+2;++i)
		{
			z[i]=new double[gpd_size];
			memset(z[i],0,sizeof(double)*gpd_size);
		}
	}

	InitAlpha();
	InitLambda();
}

// initialization of parameter alpha
void GPD::InitAlpha()
{
	if(initAlpha==NULL)
		initAlpha=new double[gpd_size];

	//	srand(time(NULL));
	for(int j=0;j<gpd_size;++j)
	{
		//		alpha[j]=(double)rand()/(double)RAND_MAX/(1+VERY_SMALL_NUMBER)/(double)gpd_size+VERY_SMALL_NUMBER;
		alpha[j]=(double)1/(double)gpd_size;
		initAlpha[j]=alpha[j];
	}
}

void GPD::InitAlpha(double* vector)
{
	alpha=vector;
}

// initialization of parameter lambda
//
void GPD::InitLambda()
{
	if(!isSimData)
	{
		// For the real data, since 0 & 1 were inserted to fit GPD,
		// 0 & 1 will be assigned to lambda, which would lead to an error.
		// Therefore, 0 & 1 must be removed.
		if(x[0]==0)
		{
			x.erase(x.begin());
			l.erase(l.begin());
		}
		if(x[0]==1)
		{
			x.erase(x.begin());
			l.erase(l.begin());
		}
	}

	// find the valid maximum occurrence of unique k-mer,
	// which will be the upper limit of lambda.
	// For the unique k-mer occurrences x[i], the numbers
	// close to upper- or lower- bound may not be valid
	// because they may be not statistically significant.
	int lambda_upper,lambda_bottom;
	int tailsum=0;
	for(int i=uniq_kmer_size-1;i>=0;--i)
	{
		tailsum+=l[i];
		if((double)tailsum/all_kmer_size>=valid_kmer_percent)
		{
			lambda_upper=x[i];
			break;
		}
	}
	int headsum=0;
	for(int i=0;i<uniq_kmer_size;++i)
	{
		headsum+=l[i];
		if((double)headsum/all_kmer_size>=valid_kmer_percent)
		{
			lambda_bottom=x[i];
			break;
		}
	}

	// Initialize lambda_j1 based on the unique k-mer occurrences,
	// while lambda_j2 is assigned a random value less than 1/gpd_size.
	double incr=(double)(lambda_upper-lambda_bottom)/gpd_size;
	for(int j=0;j<gpd_size;++j)
	{
		lambda[j][0]=lambda_bottom+(double)j*incr;
		lambda[j][1]=(double)rand()/(double)RAND_MAX/(double)(gpd_size+1)+VERY_SMALL_NUMBER;
		initLambda[j][0]=lambda[j][0];
		initLambda[j][1]=lambda[j][1];
	}
}

void GPD::InitLambda(double** matrix)
{
	lambda=matrix;
}

// EMrun - execute EM algorithm
// Inputs:
//		None.
// Outputs:
//		em_ret - 1 succeed; 0 fail
int GPD::EMrun()
{
	bool isConverge=false;
	int iter_cnt=0;
	int em_ret=1;

	// allocate memory
	double** oldLambda=new double*[gpd_size];
	for(int i=0;i<gpd_size;++i)
	{
		oldLambda[i]=new double[LAMBDA_SIZE];
	}
	double* oldAlpha=new double[gpd_size];

	// optimize the likelihood
	while(!isConverge && em_ret && iter_cnt<EM_MAX_ITER)
	{
		// Calculate z for this iteration
		for(int i=0;i<uniq_kmer_size;++i)
			for(int j=0;j<gpd_size;++j)
				z[i][j]=Calc_z(i,j);

		if(!isSimData)	// insert 0 & 1 for real data
			TuneGPD();

		// update alpha & lambda
		for(int j=0;j<gpd_size;++j)
		{
			for(int k=0;k<LAMBDA_SIZE;++k)
			{
				oldLambda[j][k]=lambda[j][k]; // record the value of lambda
			}
			oldAlpha[j]=alpha[j];

			// update alpha
			double sum_zi=0.0;
			for(int i=0;i<uniq_kmer_size;++i)
			{
				sum_zi+=l[i]*z[i][j];
			}
			alpha[j]=sum_zi/(double)all_kmer_size;

			// update lambda
			int ret=NewtonMethod(j);
			if(ret==-1 || !CheckParam(j))
			{
				// If cannot find a valid solution or GPD parameters
				// become out of bound, break the loop and return to main.
				em_ret=0;
				break;
			} // end of if
		} // end of for

		// check for convergence
		if(em_ret)
		{
			cout<<"EMrun() iter "<<iter_cnt<<":  lambda="<<endl;	// TEST ONLY
			ShowLambda();	// TEST ONLY
			//
			// converge criteria:
			//	sqrt((lambda_j1-old_lambda_j1)^2 + (lambda_j2-old_lambda_j2)^2 +
			// 			(alpha1-old_alpha1)^2 + (alpha2-old_alpha2)^2) /
			//	sqrt(old_lambda_j1^2 + old_lambda_j2^2 + old_alpha1^2 + old_alpha2^2)
			//		< em_thresh(such as 0.01)
			//
			double sum_diff=0.0;
			double sum_old=0.0;
			for(int j=0;j<gpd_size;++j)
			{
				for(int k=0;k<LAMBDA_SIZE;++k)
				{
					sum_diff=sum_diff+pow(lambda[j][k]-oldLambda[j][k],2);
					sum_old=sum_old+pow(lambda[j][k],2);
				}
				sum_diff=sum_diff+pow(alpha[j]-oldAlpha[j],2);
				sum_old=sum_old+pow(alpha[j],2);
			} // end of for
			double rate=sqrt(sum_diff)/sqrt(sum_old);
			if(rate<em_thresh)
			{
				isConverge=true;
			} else {
				isConverge=false;
			}
		} else {
			break;
		}
		iter_cnt++;
	} // end of while

	for(int i=0;i<gpd_size;++i)
	{
		delete [] oldLambda[i];
	}
	delete [] oldLambda;
	delete [] oldAlpha;

	return em_ret;
}

// insert 0 & 1 as uniq k-mers so that kmer-occurs will fit GPD better
void GPD::TuneGPD()
{
	if(isSimData)
		return;

	Classify_kmer();

	if(x[0]!=0)
	{
		x.insert(x.begin(),1);
		x.insert(x.begin(),0);
		l.insert(l.begin(),0);
		l.insert(l.begin(),0);
		uniq_kmer_size=x.size(); // enlarge uniq_kmer_size
	}

	int oldl0=l[0];
	int oldl1=l[1];
	l[0]=0;
	l[1]=0;
	for(int j=0;j<gpd_size;++j)
	{
//		double tail_0=gsl_cdf_poisson_P(0,lambda[j][0]);
//		double tail_1=gsl_cdf_poisson_P(1,lambda[j][0]);
		double tail_0=cdf(0,j);
		double tail_1=cdf(1,j);
		double new_total_kmer=kmerbin[j].kmer_size;
		for(int k=0;k<kmerbin[j].bin_size;++k)
		{
			if(kmerbin[j].kmers[k]==x[0])
				new_total_kmer-=oldl0;
			if(kmerbin[j].kmers[k]==x[1])
				new_total_kmer-=oldl1;
		}
		new_total_kmer /= (1-tail_1);
		l[0]+=(int)new_total_kmer*tail_0;
		l[1]+=(int)new_total_kmer*(tail_1-tail_0);
	}
}

// Calculate the Generalized Probability
// inputs:
//		i -- number of reads
//		j -- index of y, number of GPDs
double GPD::CalcProb(unsigned int i,unsigned int j)
{
	double ln_xi=0.0;
	for(int k=1;k<=x[i];++k)
		ln_xi += log(k);

	double prob=exp(log(lambda[j][0])+(x[i]-1)*log(lambda[j][0]+x[i]*lambda[j][1])-(lambda[j][0]+x[i]*lambda[j][1])-ln_xi);

	if(isnan(prob)||isinf(prob))
	{
		cout<<"ERROR: data out of range in GPD::CalcProb()."<<endl;
		prob=VERY_SMALL_NUMBER;
	}

	return prob;
}

// calculate probability by specifying x
double GPD::prob(int x,unsigned int j)
{
	double ln_x=0;
	for(int k=1;k<=x;++k)
		ln_x += log(k);

	double prob=exp(log(lambda[j][0])+(x-1)*log(lambda[j][0]+x*lambda[j][1])-(lambda[j][0]+x*lambda[j][1])-ln_x);
	if(isnan(prob)||isinf(prob))
	{
		cout<<"ERROR: data out of range in GPD::CalcProb()."<<endl;
		prob=VERY_SMALL_NUMBER;
	}

	return prob;
}

// cumulative distribution function
double GPD::cdf(unsigned int x,unsigned int j)
{
	double c=0.0;
	for(unsigned int i=0;i<=x;++i)
		c+=prob(i,j);

	return c;
}

// inputs:
//		i -- index of k-mers
//		j -- index of GPDs
double GPD::Calc_z(unsigned int i,unsigned int j)
{
	double denominator=0;
	for(int k=0;k<gpd_size;++k)
	{
		denominator += alpha[k]*CalcProb(i,k);
	}

	double z_tmp=alpha[j]*CalcProb(i,j)/denominator;
	if(isnan(z_tmp)||isinf(z_tmp))
	{
		//cout<<"ERROR: data out of range in GPD::Calc_z()."<<endl;
		z_tmp=VERY_SMALL_NUMBER;
	}

	return z_tmp;
}

// Group each k-mer into a cluster
void GPD::Classify_kmer()
{
	if(!isSimData)
	{
		// For the real data, since 0 & 1 were inserted to fit GPD,
		// 0 & 1 will be assigned to lambda, which would lead to an error.
		// Therefore, 0 & 1 must be removed.
		if(x[0]==0)
		{
			x.erase(x.begin());
			l.erase(l.begin());
		}
		if(x[0]==1)
		{
			x.erase(x.begin());
			l.erase(l.begin());
		}
	}

	binOfKmer.clear();
	for(int i=0;i<uniq_kmer_size;++i)
	{
		double max=Calc_z(i,0);
		binOfKmer.insert(pair<int,int>(x[i],0));
		double prob_z;
		for(int j=1;j<gpd_size;++j)
		{
			prob_z=Calc_z(i,j);
			if(max<prob_z)
			{
				max=prob_z;
				binOfKmer[x[i]]=j;
			}
		}
	} // end of for_i

	kmerbin.clear();
	// process classified k-mers
	for(int j=0;j<gpd_size;++j)
	{
		KmerBin node;
		node.bin=j;
		node.bin_size=0;
		node.kmer_size=0;
		for(int i=0;i<uniq_kmer_size;++i)
		{
			if(binOfKmer[x[i]]==j)
			{
				node.kmers.insert(node.kmers.end(),x[i]);
				node.probs.insert(node.probs.end(),CalcProb(i,j));
				node.kmer_size+=l[i];
				node.bin_size++;
			}
		}
		kmerbin.insert(kmerbin.end(),node);
	}
} // end of function Classify()

// For a k-mer, find out which bin is most probably can be classified into.
int GPD::GetMaxBin(int* a, int len)
{
	int max_bin=0;
	int max_kmer=a[0];
	for(int i=1;i<len;++i)
	{
		if(max_kmer<a[i])
		{
			max_kmer=a[i];
			max_bin=i;
		}
	}
	return max_bin;
}

// find out which bin contains the most k-mers of a read
int GPD::GetMaxBin(double* a, int len)
{
	int max_bin=0;
	double max_prob=a[0];
	for(int i=1;i<len;++i)
	{
		if(max_prob<a[i])
		{
			max_prob=a[i];
			max_bin=i;
		}
	}
	return max_bin;
}

// binning and write classified reads into logs
void GPD::WriteResult_read(char* inputFile, char* resultFile)
{
	Classify_kmer();

	fstream fs;
	char* file_name=new char[FILE_NAME_LEN];	// temp file name
	int* kmersofbin=new int[gpd_size];	// number of kmers in a read for each bin
	double* readofbin=new double[gpd_size];	// probability that a read belongs to each bin

	// calculate and record index
	double chi=CalcCHI();
	double swc=CalcSWC();

	memset(file_name,'\0',sizeof(char)*FILE_NAME_LEN);
	if(recordKmerLen)
	{
		//		sprintf(file_name,"%s_%d.score",resultFile,kmer_len);
		sprintf(file_name,"%s.score",resultFile);
	} else {
		sprintf(file_name,"%s.score",resultFile);
	}
	fs.open(file_name, fstream::out);
	if(fs.is_open() == false)
	{
		cout<<"ERROR: GPD::WriteResult_read: failed to open file "<<resultFile<<"!\n";
		return;
	}
	fs<<"Number of bins: "<<gpd_size<<endl;
	fs<<"Likelihood: "<<Likelihood()<<endl;
	fs<<"BIC: "<<BIC(LogLikelihood())<<endl;
	fs<<"CHI(Calinski-Harabasz Index): "<<chi<<endl;
	fs<<"SWC(Silhouette Width Criterion): "<<swc<<endl<<endl;
	fs<<"Final values of lambdas:"<<endl;
	for(int j=0;j<gpd_size;++j)
	{
		for(int k=0;k<LAMBDA_SIZE;++k)
			fs<<lambda[j][k]<<" ";
		fs<<endl;
	}
	fs<<endl;
	fs<<"Final values of alphas:"<<endl;
	for(int j=0;j<gpd_size;++j)
	{
		fs<<alpha[j]<<" ";
	}
	fs<<endl;
	fs<<"Convergence criteria(z=E(y)):\n EM rate < "<<em_thresh<<"; Newton rate < "<<NEWTON_EPS<<endl<<endl;
	fs.close();

	//return; //// TEST ONLY ////

	gzFile fp;
	kseq_t *seq; // genome sequence
	int l=0;

	fp = gzopen(inputFile, "r"); // STEP 2: open the file handler
	seq = kseq_init(fp); // STEP 3: initialize seq
	while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
		memset(kmersofbin,0,sizeof(int)*gpd_size);
		//vector<int> kmersInRead;
		// read k-mers
		for(unsigned int i=0;i<strlen(seq->seq.s)-kmer_len;++i)
		{
			string k_mer(seq->seq.s+i,kmer_len);
			map<string,int>::iterator it;
			it=kmer_occur.find(k_mer);
			if((*it).second>1 && (*it).second<=x[uniq_kmer_size-1])	// the k-mers occurs only once was ignored
			{
				++kmersofbin[binOfKmer[(*it).second]];
				//kmersInRead.insert(kmersInRead.end(),(*it).second);
			}
		} // end of for

		int bin_no=GetMaxBin(kmersofbin,gpd_size);
		/*for(int j=0;j<gpd_size;++j)
		{
			readofbin[j]=CalcReadProb(kmersInRead,j);
		}
		int bin_no=GetMaxBin(readofbin,gpd_size);*/

		// write this bin into corresponding file
		memset(file_name,'\0',sizeof(char)*FILE_NAME_LEN);
		if(recordKmerLen)
		{
			//sprintf(file_name,"%s_kmer%d.%d",resultFile,kmer_len,bin_no);
			sprintf(file_name,"%s.%d",resultFile,bin_no);
		} else {
			sprintf(file_name,"%s.%d",resultFile,bin_no);
		}

		fs.open(file_name, fstream::out|fstream::app);
		if(fs.is_open() == false)
		{
			cout<<"ERROR: GPD::WriteResult_read: failed to open file "<<resultFile<<"!\n";
			return;
		}
		fs<<">"<<seq->name.s<<" ";
		if(seq->comment.s!=NULL && strlen(seq->comment.s)>0)
			fs<<seq->comment.s<<endl;
		else
			fs<<endl;
			fs<<seq->seq.s<<endl;
		fs.close();
	}
	kseq_destroy(seq); // STEP 5: destroy seq
	gzclose(fp); // STEP 6: close the file handler

	delete [] kmersofbin;
	delete [] file_name;
	delete [] readofbin;
}

// calculate the probability that a read belongs to a bin
double GPD::CalcReadProb(vector<int> &kmers, int j)
{
	double p;
	double numer=1;
	double denorm=0.0;
	for(unsigned int i=0;i<kmers.size();++i)
	{
		numer *= Calc_z(GetIndexofx(kmers[i]),j);
	}
	for(int k=0;k<gpd_size;++k)
	{
		double prod=1;
		for(unsigned int i=0;i<kmers.size();++i)
		{
			prod *= Calc_z(GetIndexofx(kmers[i]),k);
		}
		denorm+=prod;
	}

	p=numer/denorm;

	return p;
}

// Get the index of kmer in x
int GPD::GetIndexofx(int occur)
{
	int ret=-1;
	for(int i=0;i<uniq_kmer_size;++i)
	{
		if(occur==x[i])
		{
			ret=i;
			break;
		}
	}
	return ret;
}

// This function can only be called after Classify_().
void GPD::WriteResult_kmer(char* resultFile)
{
	Classify_kmer();

	if(binOfKmer.empty())
	{
		cout<<"ERROR: reads have not been classified!"<<endl;
		exit(1);
	}

	fstream fs;
	fs.open(resultFile, fstream::out|fstream::app);
	if(fs.is_open() == false)
	{
		cout<<"ERROR: GPD::WriteResult_kmer: failed to open file "<<resultFile<<"!\n";
		return;
	}

	double chi=CalcCHI();
	double swc=CalcSWC();
	fs<<"\n================================================\n";
	fs<<"Calinski-Harabasz Index: "<<chi<<endl;
	fs<<"Silhouette Width Criterion: "<<swc<<endl;

	// Record initial value of lambda, so that we can estimate the range of workable lambda.
	fs<<"Initial values of lambdas are:"<<endl;
	for(int j=0;j<gpd_size;++j)
	{
		for(int k=0;k<LAMBDA_SIZE;++k)
			fs<<initLambda[j][k]<<" ";
		fs<<endl;
	}
	fs<<endl;

	fs<<"Final values of lambdas are:"<<endl;
	for(int j=0;j<gpd_size;++j)
	{
		for(int k=0;k<LAMBDA_SIZE;++k)
			fs<<lambda[j][k]<<" ";
		fs<<endl;
	}
	fs<<endl;

	fs<<"Initial values of alphas are:"<<endl;
	for(int j=0;j<gpd_size;++j)
	{
		fs<<initAlpha[j]<<" ";
	}
	fs<<endl;

	fs<<"Final values of alphas are:"<<endl;
	for(int j=0;j<gpd_size;++j)
	{
		fs<<alpha[j]<<" ";
	}
	fs<<endl;

	fs<<"Convergence criteria(z=E(y)):\n EM rate < "<<CONVERGE_THRESHHOLD<<"; Newton rate < "<<NEWTON_EPS<<endl<<endl;

	int* total_bin_size=new int[gpd_size];
	memset(total_bin_size,0,sizeof(int)*gpd_size);
	for(int j=0;j<gpd_size;++j)
	{
		for(int i=0;i<uniq_kmer_size;++i)
		{
			if(binOfKmer[x[i]]==j)
			{
				total_bin_size[j]+=l[i];
			}
		} // end of for_i
	}

	//	double total_diff=0.0;
	// Log the result of GPD.
	for(int j=0;j<gpd_size;++j)
	{
		fs<<"Generalized Poisson Distribution "<<j<<"(total "<<total_bin_size[j]<<"):"<<endl;
		for(int i=0;i<uniq_kmer_size;++i)
		{
			if(binOfKmer[x[i]]==j)
			{
				//				for(int k=0;k<l[i];++k)
				fs<<x[i]<<" ";
			}
		} // end of for_i
		fs<<endl;
	} // end of for_j

	fs.close();

	delete [] total_bin_size;
}
// Print lambdas. TEST ONLY
void GPD::ShowLambda()
{
	for(int j=0;j<gpd_size;++j)
	{
		for(int k=0;k<LAMBDA_SIZE;++k)
			cout<<lambda[j][k]<<" ";
		cout<<endl;
	}
}

double GPD::Calc_w(int j)
{
	double denorm=0.0;
	double numer=0.0;

	for(int i=0;i<uniq_kmer_size;++i)
	{
		denorm += z[i][j]*l[i];
		numer += z[i][j]*x[i]*l[i];
	}
	if(isnan(denorm)||isnan(numer) || isinf(denorm)||isinf(numer))
	{
		cout<<"ERROR: data out of range in GPD::Calc_w()."<<endl;
	}
	return numer/denorm;
}

double GPD::f_lambda_j2(int j)
{
	double w=Calc_w(j);
	double sum=0.0;

	for(int i=0;i<uniq_kmer_size;++i)
	{
		sum = sum + l[i]*z[i][j]*(x[i]*(x[i]-1)/(w+(x[i]-w)*lambda[j][1])-x[i]);
	}

	if(isnan(sum)||isinf(sum))
	{
		cout<<"ERROR: data out of range in GPD::f_lambda_j2()."<<endl;
	}

	return sum;
}

double GPD::df_lambda_j2(int j)
{
	double w=Calc_w(j);
	double sum=0.0;

	for(int i=0;i<uniq_kmer_size;++i)
	{
		sum = sum - (l[i]*z[i][j]*x[i]*(x[i]-1)*(x[i]-w))/pow(w+(x[i]-w)*lambda[j][1],2);
	}

	if(isnan(sum)||isinf(sum))
	{
		cout<<"ERROR: data out of range in GPD::df_lambda_j2()."<<endl;
	}

	return sum;
}

// Calculate lambda_j1 & lambda_j2
int GPD::NewtonMethod(int j)
{
	double f,df;
	double oldLambda_j2;
	int cnt=0;
	int ret=0;

	// Calculate lambda_j1
	while(cnt<NEWTON_MAX_ITER)
	{
		oldLambda_j2=lambda[j][1];

		f=f_lambda_j2(j);
		df=df_lambda_j2(j);
		lambda[j][1]-=f/df;

		double rate=abs((lambda[j][1]-oldLambda_j2)/lambda[j][1]);
		if(isnan(rate) || isinf(rate) || !CheckParam(j))
		{
			ret=-1;
			break;
		} else if(rate<NEWTON_EPS)
		{
			break;
		}
		cnt++;
	}

	// Calculate lambda_j2
	if(!ret)
		lambda[j][0]=Calc_w(j)*(1-lambda[j][1]);

	return ret;
}

// Check the validity of alphas and lambdas
int GPD::CheckParam(int j)
{
	int ret=1;
	if(lambda[j][1]<0 || lambda[j][1]>1 || alpha[j]<0 || alpha[j]>1)
	{
		ret=0;
	}

	return ret;
}

// calculate Calinski-Harabasz Index
double GPD::CalcCHI()
{
	// calculate between-cluster and within-cluster sum of squares.
	double chi=0.0;
	double b=0.0;	// between-cluster sum
	double w=0.0;	// with-cluster sum

	double* cent=new double[gpd_size];
	memset(cent,0,sizeof(double)*gpd_size);

	// within-cluster sum of squares
	for(int j=0;j<gpd_size;++j)
	{
		cent[j]=0.0; // cluster center
		for(int i=0;i<kmerbin[j].bin_size;++i)
		{
			cent[j]+=kmerbin[j].probs[i];
		}
		if(kmerbin[j].bin_size!=0)
			cent[j]/=kmerbin[j].bin_size;
		else
			cent[j]=0;

		for(int i=0;i<kmerbin[j].bin_size;++i)
		{
			w+=pow(kmerbin[j].probs[i]-cent[j],2);
		}
	}

	// calc the centroid of the entire dataset
	int datasize=0;
	double centroid=0.0;
	for(int j=0;j<gpd_size;++j)
	{
		datasize+=kmerbin[j].bin_size;
		for(int i=0;i<kmerbin[j].bin_size;++i)
		{
			centroid+=kmerbin[j].probs[i];
		}
	}
	centroid /= datasize;

	// between-cluster sum of squares
	for(int j=0;j<gpd_size;++j)
	{
		b+=kmerbin[j].bin_size*pow(cent[j]-centroid,2);
	} // end of for j

	delete [] cent;

	chi=(b/(gpd_size-1))/(w/(datasize-gpd_size));
	return chi;
}

// calculate Silhouette Width Criterion
double GPD::CalcSWC()
{
	// a[j][i] is the average distance of this object to all other objects
	// in the same cluster.
	double* a=new double[uniq_kmer_size];
	// b[j][i] is the minimum distance of this object to all objects in another cluster.
	double* b=new double[uniq_kmer_size];
	double* d=new double[gpd_size-1];
	double* s=new double[uniq_kmer_size];
	memset(a,0,sizeof(double)*uniq_kmer_size);
	memset(b,0,sizeof(double)*uniq_kmer_size);

	for(int i=0;i<uniq_kmer_size;++i)
	{
		memset(d,0,sizeof(double)*(gpd_size-1));
		int idx=0;
		for(int j=0;j<gpd_size;++j)
		{
			if(binOfKmer[x[i]]==kmerbin[j].bin)
			{
				// calculate intra-cluster distance
				if(kmerbin[j].bin_size==1)
				{
					a[i]=0;
				} else {
					for(int k=0;k<kmerbin[j].bin_size;++k)
					{
						a[i]+=CalcJSD(x[i],kmerbin[j].kmers[k]);
					}
					a[i]/=(double)(kmerbin[j].bin_size-1);
				}
			} else {
				// calculate inter-cluster distance
				for(int len=0;len<kmerbin[j].bin_size;++len)
				{
					d[idx]+=CalcJSD(x[i],kmerbin[j].kmers[len]);
				}
				d[idx]/=(double)kmerbin[j].bin_size;
				idx++;
			}
		} // end of for j
		double min=d[0];
		for(int p=1;p<idx;++p)
		{
			if(min>d[p])
				min=d[p];
		}
		b[i]=min;

		s[i]=(b[i]-a[i])/(b[i]>a[i]?b[i]:a[i]);
		/*		if(isnan(s[i])||isinf(s[i]))
		{
			cout<<"ERROR: invalid score of SWC."<<endl;
			cout<<"a[i]="<<a[i]<<endl;
			cout<<"b[i]="<<b[i]<<endl;
		}*/
	} // end of for i

	double swc=0.0;
	for(int i=0;i<uniq_kmer_size;++i)
		swc+=s[i];
	swc/=(double)uniq_kmer_size;

	// release allocated memory
	delete [] a;
	delete [] b;
	delete [] d;
	delete [] s;

	return swc;
}

// distance between points (x1,p1) and (x2,p2)
double GPD::dist(int x1,double p1,int x2,double p2)
{
	return sqrt(pow((double)x1-(double)x2,2)+pow(p1-p2,2));
}

double GPD::dist(double x,double y)
{
	return sqrt(pow(x-y,2));
}

// calculate Kullback-Leibler Divergence(KLD)
double GPD::CalcKLD(double* x,double* y,int len)
{
	double kld=0.0;
	for(int i=0;i<len;++i)
	{
		kld+=x[i]*log(x[i]/y[i]);
	}
	return kld;
}

// calculate Kullback-Leibler Divergence(KLD)
double GPD::CalcKLD(double x,double y)
{
	if(x==0)
		x+=0.0000001;
	if(y==0)
		y+=0.0000001;

	double ret=x*log(x/y);
	if(isnan(ret)||isinf(ret))
	{
		cout<<"ERROR: invalid KLD!"<<endl;
	}

	return ret;
}

// calculate distance
// Jensen-Shannon Divergence(JSD) and Kullback-Leibler Divergence(KLD) and used.
// inputs:
//		x,y - arrays containing the probabilities of two distributions.
double GPD::CalcJSD(double* x,double* y,int len)
{
	double jsd=0.0;
	double* m=new double[len];
	for(int i=0;i<len;++i)
		m[i]=(x[i]+y[i])/2;

	jsd=(CalcKLD(x,m,len)+CalcKLD(y,m,len))/2;

	delete [] m;
	return sqrt(jsd);
}

double GPD::CalcJSD(double x,double y)
{
	double m=(x+y)/2;
	double jsd=(CalcKLD(x,m)+CalcKLD(y,m))/2;
	if(isnan(jsd)||isinf(jsd))
	{
		cout<<"ERROR: invalid JSD!"<<endl;
	}

	return jsd;
}

// cummulative hypergeometric distribution function
// The probability distribution for hypergeometric random variates is,
//
//         p(k) =  C(n_1, k) C(n_2, t - k) / C(n_1 + n_2, t)
// where C(a,b) = a!/(b!(a-b)!) and t <= n_1 + n_2. The domain of k is
// max(0,t-n_2), ..., min(t,n_1).
//
// If a population contains n_1 elements of “type 1” and n_2 elements of “type 2”
// then the hypergeometric distribution gives the probability of obtaining k
// elements of “type 1” in t samples from the population without replacement.
double GPD::phyper(unsigned int k, unsigned int n1, unsigned int n2, unsigned int t)
{
	double ret=gsl_cdf_hypergeometric_P(k,n1,n2,t);;

	return ret;
}

void GPD::WriteScore(char* resultFile)
{
	fstream fs;
	char* file_name=new char[FILE_NAME_LEN];	// temp file name

	// calculate and record index
	double chi=CalcCHI();
	double swc=CalcSWC();
	double ii=IndexI();

	memset(file_name,'\0',sizeof(char)*FILE_NAME_LEN);
	sprintf(file_name,"%s.score.%d",resultFile,gpd_size);
	fs.open(file_name, fstream::out);
	if(fs.is_open() == false)
	{
		cout<<"ERROR: failed to open file "<<resultFile<<"!\n";
		return;
	}
	fs<<"Number of bins: "<<gpd_size<<endl;
	fs<<"CHI(Calinski-Harabasz Index): "<<chi<<endl;
	fs<<"SWC(Silhouette Width Criterion): "<<swc<<endl;
	fs<<"Index I: "<<ii<<endl;
	fs.close();

	delete [] file_name;
}

double GPD::Likelihood()
{
	double lh=1;
	double sum=0;

	for(int j=0;j<gpd_size;++j)
	{
		for(int i=0;i<uniq_kmer_size;++i)
		{
			sum += z[i][j]*alpha[j]*CalcProb(i,j);
		}
		lh *= sum;
	}

	return lh;
}

// Calculate the likelihood of L(X,Y,L|Theta)
double GPD::LogLikelihood()
{
	double llh=0;
	// approximate the likelihood by log likelihood
	for(int i=0;i<uniq_kmer_size;++i)
	{
		double ln_x=0;
		for(int k=1;k<=x[i];++k)
			ln_x += log(k);

		for(int j=0;j<gpd_size;++j)
		{
			llh += z[i][j]*(log(alpha[j])+log(lambda[j][0])+(x[i]-1)*
					log(lambda[j][0]+x[i]*lambda[j][1])-(lambda[j][0]+x[i]*lambda[j][1])-ln_x);
		}
	}

	return llh;
}

// Calculate the Bayesian Information Criterion
double GPD::BIC(double lh)
{
	return -2*lh+(gpd_size*3)*log(uniq_kmer_size);
}

double GPD::IndexI()
{
	// find centers
	double* cent=new double[gpd_size];
	for(int j=0;j<gpd_size;++j)
	{
		cent[j]=0.0; // cluster center
		for(int i=0;i<kmerbin[j].bin_size;++i)
		{
			cent[j]+=kmerbin[j].probs[i];
		}
		if(kmerbin[j].bin_size!=0)
		{
			cent[j]/=kmerbin[j].bin_size;
		} else {
			cent[j]=0;
		}
	}

	// Calculate E1
	double E1=0;
	for(int i=0;i<kmerbin[0].bin_size;++i)
	{
		E1+=sqrt(pow(kmerbin[0].probs[i]-cent[0],2));
	}

	// Calculate EK
	double EK=0;
	for(int j=0;j<gpd_size;++j)
	{
		double sum=0;
		for(int i=0;i<kmerbin[j].bin_size;++i)
		{
			sum+=sqrt(pow(kmerbin[j].probs[i]-cent[j],2));
		}
		EK+=sum;
	}

	// Calcualte DK
	double DK=-1;
	for(int j=0;j<gpd_size-1;++j)
	{
		for(int i=j+1;i<gpd_size;++i)
		{
			double d=sqrt(pow(cent[i]-cent[j],2));
			if(DK<d)
				DK=d;
		}
	}

	delete [] cent;

	double score=pow(((double)1/(double)gpd_size*E1/EK*DK),2);
	return score;
}
