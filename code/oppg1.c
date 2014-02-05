#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "common.h"

double sum(int N);
double sumSlow(int N);
double sumShared(int N);

int main(int argc,char** argv){
	double S = (M_PI*M_PI)/6;
	int iterations = 14, i;
	int N[iterations];
	
	
	N[0] = 8;
	for(i=1; i<iterations; ++i)
	{
		N[i] = 2*N[i-1];
	}
	double* Sn = (double*)malloc(iterations*sizeof(double));
	double* SnSlow = (double*)malloc(iterations*sizeof(double));
	
	for(i=0; i<iterations; ++i)
	{
		Sn[i] = sumShared(N[i]);
		SnSlow[i] = sum(N[i]);
	}
	for(i=0; i<iterations; ++i)
	{
		printf("%e %e\n",S-SnSlow[i], Sn[i]-SnSlow[i]);
	}
	free(Sn); free(SnSlow);
	return 0;
}

// Makes and sums up a vector. Sums in the numerically worst possible way.
double
sum(int N)
{
	int i;
	double Sn=0;
	double* v = (double*)malloc(N*sizeof(double));
	for(i=1;i<=N;++i)
	{
		v[i-1] = 1.0/(i*i);
	}
	for(i=0; i<N; ++i)
	{
		Sn += v[i];
	}
	return Sn;
}

// Slower edition of sum. Sums by always adding the two smallest numbers. Numarically optimal.
double
sumSlow(int N)
{
	int i, j;
	double* v = (double*)malloc(N*sizeof(double));
	double tmp;
	for(i=1;i<=N;++i)
	{
		v[i-1] = 1.0/(i*i);
	}
	for(i=N-1; i>=1; --i)
	{
		v[i-1] += v[i];
		j=i-1;
		while(v[j]>v[j-1] && j>0)
		{
			tmp = v[j-1];
			v[j-1] = v[j];
			v[j] = tmp;
			--j;
		}
	}
	tmp = v[0];
	free(v);
	return tmp;
}
// Shared memory parallelisation.
double
sumShared(int N)
{
	int i;
	double Sn=0;
	double* v = (double*)malloc(N*sizeof(double));

#pragma omp parallel for schedule(static)
	for(i=1; i<=N; ++i)
	{
		v[i-1] = 1.0/(i*i);
	}

#pragma omp parallel for schedule(static) reduction(+:Sn)
	for(i=0; i<N; ++i)
	{
		Sn += v[i];
	}
	return Sn;
}

