#include <stdlib.h>
#include "sum.h"

// Makes and sums up a vector. Sums from lowest to highest element.
double
sum(const long long N)
{
	double* v = initVec_seq(N);
	double Sn = summingSeq(v, N);
	free(v);
	return Sn;
}

// Makes and sums up a vector. Sums from highest to lowest element.
double
sumF(const long long N)
{
	double* v = initVec_seq(N);
	double Sn = summingSeqF(v, N);
	free(v);
	return Sn;
}

// Slower edition of sum. Sums by always adding the two smallest numbers.
// Numarically more stable.
double
sumSlow(const long long N)
{
	double* v = initVec_seq(N);
	double Sn = summingSlow(v, N);
	free(v);
	return Sn;
}
// Shared memory parallelisation.
double
sumShared(const long long N)
{
	double* v = initVec_shared(N);
	double Sn = summingShared(v, N);
	free(v);
	return Sn;
}

// Distributed memory parallelisation.
double
sumDist(const long long N, int* rank)
{
	double Sn=0;
	double Sn_temp;
	int size;
	long long N_part;
	double* v;

	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);

	//The following line is OK only as long as N = 2^k1 and P=2^k2 and k1>k2.
	N_part = N/size;
	v = initVec_dist(N_part, rank, size);
	Sn_temp = summingSeq(v,N_part);
	free(v);
	
	MPI_Reduce(&Sn_temp,&Sn,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);	
	if(*rank == 0)
		return Sn;
	else
		return 0.0;
}

// Parallelisation using both MPI and OpenMP.
double
sumHybrid(const long long N, int* rank)
{
	double Sn=0;
	double Sn_temp;
	int size;
	long long N_part;
	double* v;

	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);

	//The following line is OK only as long as N = 2^k1 and P=2^k2 and k1>k2.
	N_part = N/size;
	v = initVec_hybrid(N_part, rank, size);
	Sn_temp = summingShared(v,N_part);
	free(v);
	
	MPI_Reduce(&Sn_temp,&Sn,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);	
	if(*rank == 0)
		return Sn;
	else
		return 0.0;
}

