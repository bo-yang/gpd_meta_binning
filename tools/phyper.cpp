#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>
#include <gsl/gsl_cdf.h>

using namespace std;

// 
// Build options:
//
// g++ -o phyper phyper.cpp -lgsl -lgslcblas -lm
//

int main(int argc, char* argv[])
{
	unsigned long int k;
	unsigned long int n1;
	unsigned long int n2;
	unsigned long int t;
	
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-k") == 0)
		{
			i++;
			k = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-n1") == 0)
		{
			i++;
			n1 = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-n2") == 0)
		{
			i++;
			n2 = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-t") == 0)
		{
			i++;
			t = atoi(argv[i]);
		} else {
			cout<<"ERROR: invalid parameters "<<argv[i]<<endl;
			cout<<"Usage: phyper -k <k> -n1 <n1> -n2 <n2> -t <t>."<<endl;
			return 1;
		}
	}

	// Cumulative Poisson Distribution
	double p=gsl_cdf_hypergeometric_P(k,n1,n2,t);
	//double q=gsl_cdf_hypergeometric_Q(k,n1,n2,t);
	cout<<setprecision(12)<<p<<endl;
	//cout<<setprecision(12)<<q<<endl;

	return 0;
}
