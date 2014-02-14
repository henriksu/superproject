#include "sum.h"

// Sums the N elements of the array v from right to left.
double summingSeq(double* v, const long long N)
{
	long long i;
	double Sn = 0;
	for(i=N-1; i>=0; --i)
	{
		Sn += v[i];
	}
	return Sn;
}

// Sums the N elements of the array v from left to right.
double summingSeqF(double* v, const long long N)
{
	long long i;
	double Sn = 0;
	for(i=0; i<N; ++i)
	{
		Sn += v[i];
	}
	return Sn;
}

// Sums the N elements of the array v using OpenMP to sum different parts
// of the array. The order of summation (and thus the exact result) is
// non-deterministic.
double summingShared(double* v, const long long N)
{
	long long i;
	double Sn = 0;
	#pragma omp parallel for schedule(static) reduction(+:Sn)
	for(i=N-1; i>=0; --i)
	{
		Sn += v[i];
	}
	return Sn;
}

// Sums the N elements of the array v. Assuming v is sorted, this function will
// always sum the two smalles elements and store the result to the array.
// The approach enhanches numerical stability at the cost of sorting.
// NOTE: The array is wrecked.
double summingSlow(double* v, const long long N)
{
	long long i, j;
	double tmp;
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
	return v[0];
}

