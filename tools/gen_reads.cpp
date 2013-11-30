#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

//
// Command to build gen_reads(GSL is required):
//	g++ -Wall -o gen_reads gen_reads.cpp -lgsl -lgslcblas -lm
//

const long int GENOME_MAX_LEN=20000000;
string gen;
string gen_name;

int gen_poisson (int* posnum, int n, double mu)
{
	const gsl_rng_type * T;
	gsl_rng * r;

	/* create a generator chosen by the 
	 * environment variable GSL_RNG_TYPE */

	gsl_rng_env_setup();

  	T = gsl_rng_default;
  	r = gsl_rng_alloc (T);

  	/* generate n random variates chosen from 
	 * the poisson distribution with mean 
	 * parameter mu */

  	for (int i = 0; i < n; i++) 
    	{
      		posnum[i] = gsl_ran_poisson (r, mu);
    	}

  	gsl_rng_free (r);
  	return 0;
}

int gen_random(int* posnum, int n, int m)
{
	for(int i=0;i<n;++i)
	{
		posnum[i]=rand()%m;
	}
	return 0;
}

int load_genome(const char* genomeFile)
{
	fstream fs;
	fs.open(genomeFile, fstream::in);

	char* buf=new char[1024];
	while(!fs.eof())
	{
		fs.getline(buf,1024);
		if(buf[0]=='>')
		{
			gen_name=buf;
			gen_name.erase(0,1); // delete leading '>'
		} else {
			gen+=buf;
		}
	}

	fs.close();
	delete [] buf;

	return 0;
}

int write_reads(const char* outFile)
{
	return 0;
}

void Merge(int* input, long p, long r)
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
void Merge_sort(int* input, long p, long r)
{
	if ( p < r )
	{
		long mid = floor((p + r) / 2);
		Merge_sort(input, p, mid);
		Merge_sort(input, mid + 1, r);
		Merge(input, p, r);
	}
}

//
// Main Starts Here
//

int main(int argc, char* argv[])
{
	string input;
	string output;
	long int G = GENOME_MAX_LEN;	// length of genome
	double C = 2; 	// coverage
	unsigned int l = 100; 	// read length
	unsigned int np = 1; 	// # of poisson distributions per genome
	unsigned int nd; 	// # of sample data per Poisson distribution
	unsigned int nr; 	// # of reads
	int rd_suffix=0;	// read ID suffix, such as >r1234.1

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-i") == 0)
		{
			i++;
			input = argv[i];
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			i++;
			output = argv[i];
		}
		else if (strcmp(argv[i], "-c") == 0)
		{
			i++;
			C = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-l") == 0)
		{
			i++;
			l = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-s") == 0)
		{
			i++;
			rd_suffix = atoi(argv[i]);
		}
	/*	else if (strcmp(argv[i], "-g") == 0)
		{
			i++;
			G = atoi(argv[i]);
		}*/
		else if (strcmp(argv[i], "-np") == 0)
		{
			i++;
			np = atoi(argv[i]);
		} else {
			cout<<"ERROR: invalid parameters "<<argv[i]<<endl;
			cout<<"Usage: gen_reads -i <input-genome> -o <output> -c <coverage> -s <suffix_of_read_ID> -l <read_length>."<<endl;
			return 1;
		}
	}

	if(input.empty())
	{
		cout<<"ERROR: no genome file specified!"<<endl;
		return 1;
	}

	if(output.empty())
	{
		output="gpdreads_"+input;
	}

	load_genome(input.c_str());
	G=gen.size(); 

	// slice genome and write reads
	nr=G*C/l;
	nd=nr/np;
	srand(time(NULL));
	int* posnum=new int[nd+1]; // Start position of each read

	fstream fs;
	fs.open(output.c_str(), fstream::out|fstream::app);

	long int cnt=0;
	for(unsigned int i=0;i<np;++i)
	{
	//	unsigned int mu=rand() % (int)(G*0.8)+(int)(G*0.05);
	//	gen_poisson(posnum,nd,mu);
		gen_random(posnum,nd,G-l+1);

		Merge_sort(posnum,0,nd);
		
		for(unsigned int i=0;i<nd;++i)
		{
			fs<<">r"<<cnt<<"."<<rd_suffix<<" |"<<"LOC="<<posnum[i]<<"|"<<gen_name.c_str()<<endl;
			fs<<gen.substr(posnum[i],l)<<endl;
			cnt++;
		}
	}

	fs.close();
	delete [] posnum;
	posnum=NULL;
	
	return 0;
}
