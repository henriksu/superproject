#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "common.h"
//#include "sum.h"

double sum(const int N);
double sumSlow(const int N);
double sumShared(const int N);
double sumDist(const int N, int* rank);
double sumHybrid(const int N, int* rank);

double* initVec_seq(const int N);
double* initVec_shared(const int N);
double* initVec_dist(const int N_part, const int* rank, const int size);
double* initVec_hybrid(const int N_part, const int* rank, const int size);

double summingSeq(double* v, const int N);
double summingSlow(double* v, const int N);
double summingShared(double* v, const int N);


void basicSetup(const int iterations, int N[], double* Sn, double* SnSlow); // Fills N[] with "iterations" powers of 2, starting at 8.


int main(int argc,char** argv){
	int rank = 0;
	int i; // Generic loop variable.
	double startTime; // Storing the start time while measuring.
	double S = (M_PI*M_PI)/6; // The limit of the series.
	MPI_Init(&argc,&argv);

	// Setting up the size of the partial sums to generate. This should be altered to read something from the command line.
	int iterations = 13; // Number of different summing lengths.
	int N[iterations]; // Vector with the summetion lengths.
	double* Sn     = (double*)malloc(iterations*sizeof(double));
	double* SnSlow = (double*)malloc(iterations*sizeof(double)); // Vectors of the partial sums.
	basicSetup(iterations, N, Sn, SnSlow);

	
	for(i=0; i<iterations; ++i)
	{
		Sn[i] = sum(N[i]);//, &rank);
		SnSlow[i] = sumSlow(N[i]);
	}

	if(rank==0)
		for(i=0; i<iterations; ++i)
		{
			printf("%e %e\n",S-SnSlow[i], Sn[i]-SnSlow[i]);
		}
	free(Sn); free(SnSlow);
	
	MPI_Finalize();
	return 0;
}

void
basicSetup(const int iterations, int *N, double* Sn, double* SnSlow)
{
	int i;
	N[0] = 8;
	for(i=1; i<iterations; ++i)
	{
		N[i] = 2*N[i-1];
	}
	return;
}

// Makes and sums up a vector. Sums in the numerically worst possible way.
double
sum(const int N)
{
	double* v = initVec_seq(N);
	double Sn = summingSeq(v, N);
	free(v);
	return Sn;
}

// Slower edition of sum. Sums by always adding the two smallest numbers. Numarically optimal.
double
sumSlow(const int N)
{
	double* v = initVec_seq(N);
	double Sn = summingSlow(v, N);
	free(v);
	return Sn;
}
// Shared memory parallelisation.
double
sumShared(const int N)
{
	double* v = initVec_shared(N);
	double Sn = summingShared(v, N);
	free(v);
	return Sn;
}

double
sumDist(const int N, int* rank)
{
	double Sn=0;
	double Sn_temp;
	int size;
	int N_part;
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

double
sumHybrid(const int N, int* rank)
{
	double Sn=0;
	double Sn_temp;
	int size;
	int N_part;
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

double* initVec_seq(const int N)
{
	int i;
	double* v = (double*)malloc(N*sizeof(double));
	for(i=1; i<=N; ++i)
	{
		v[i-1] = 1.0/(i*i);
	}
	return v;
}

double* initVec_shared(const int N)
{
	int i;
	double* v = (double*)malloc(N*sizeof(double));
	#pragma omp parallel for schedule(static)
	for(i=1; i<=N; ++i)
	{
		v[i-1] = 1.0/(i*i);
	}
	return v;
}

double* initVec_dist(const int N_part, const int* rank, const int size)
{
	MPI_Status Status;
	double* v = (double*)malloc(N_part*sizeof(double));
	if(*rank == 0)
	{
		int dest,j,move;		
		for(dest=1; dest<=size; ++dest) // The last iteration is not sent to rank==size, but stays in rank==0.
		{
			move = (dest-1)*N_part+1;
			for(j=0;j<N_part;++j)
				v[j] = 1.0/((j+move)*(j+move));
			if(dest < size)
				MPI_Send(v,N_part,MPI_DOUBLE,dest,100,MPI_COMM_WORLD);
		}
	}else{
		MPI_Recv(v,N_part,MPI_DOUBLE,0,100,MPI_COMM_WORLD,&Status);
	}
	return v;
}

double* initVec_hybrid(const int N_part, const int* rank, const int size)
{
	MPI_Status Status;
	double* v =  (double*)malloc(N_part*sizeof(double));
	if(*rank == 0)
	{
		int dest,j, move = 1;
		for(dest=1; dest<=size; ++dest) // The last iteration is not sent to rank==size, but stays in rank==0.
		{
			#pragma parallel for schedule(Static)
			for(j=0;j<N_part;++j)
				v[j] = 1.0/((j+move)*(j+move));
			if(dest < size)
				MPI_Send(v,N_part,MPI_DOUBLE,dest,100,MPI_COMM_WORLD);
			move += N_part;
		}
	}else{
		MPI_Recv(v,N_part,MPI_DOUBLE,0,100,MPI_COMM_WORLD,&Status);
	}
	return v;
}


double summingSeq(double* v, const int N)
{
	int i;
	double Sn = 0;
	for(i=0; i<N; ++i)
	{
		Sn += v[i];
	}
	return Sn;
}

double summingSlow(double* v, const int N)
{
	int i, j;
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

double summingShared(double* v, const int N)
{
	int i;
	double Sn = 0;
	#pragma omp parallel for schedule(static) reduction(+:Sn)
	for(i=0; i<N; ++i)
	{
		Sn += v[i];
	}
	return Sn;
}
