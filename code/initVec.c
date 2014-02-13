#include <stdlib.h>
#include "sum.h"

// Creates an array of the N first elements of the sequanece 1/(i*i).
// Straight forward sequential implementation.
// Returns a pointer to the first element of the array.
double* initVec_seq(const long long N)
{
	long long i;
	double* v = (double*)malloc(N*sizeof(double));
	for(i=1; i<=N; ++i)
	{
		v[i-1] = 1.0/((double)i*i);
	}
	return v;
}

// Creates an array of the N first elements of the sequence 1/(i*i).
// Uses OpenMP to calculate and assign the element values in parallel.
// Returns a pointer to the first element of the array.
double* initVec_shared(const long long N)
{
	long long i;
	double* v = (double*)malloc(N*sizeof(double));
	#pragma omp parallel for schedule(static)
	for(i=1; i<=N; ++i)
	{
		v[i-1] = 1.0/((double)i*i);
	}
	return v;
}

// Create an array of consecutive elements from the sequence 1/(i*i) on each
// process. Process 0 calculates all the elements and send finished arrays
// to the other processors. Process 0 keeps the last chunk.
double* initVec_dist(const long long N_part, const int* rank, const int size)
{
	MPI_Status Status;
	double* v = (double*)malloc(N_part*sizeof(double));
	if(*rank == 0)
	{
		int dest;
		long long j, move = 1;		
		for(dest=1; dest<=size; ++dest) // The last iteration is not sent to rank==size, but stays in rank==0.
		{
			move = (dest-1)*N_part+1;
			for(j=0;j<N_part;++j)
				v[j] = 1.0/((double)(j+move)*(j+move));
			if(dest < size)
				MPI_Send(v,N_part,MPI_DOUBLE,dest,100,MPI_COMM_WORLD);
			move += N_part;
		}
	}else{
		MPI_Recv(v,N_part,MPI_DOUBLE,0,100,MPI_COMM_WORLD,&Status);
	}
	return v;
}

// Create an array of consecutive elements from the sequence 1/(i*i) on each
// process. Process 0 calculates all the elements and send finished arrays
// to the other processors. Process 0 keeps the last chunk.
// On process 0 OpenMP is used to fill each chunk in parallel.
double* initVec_hybrid(const long long N_part, const int* rank, const int size)
{
	MPI_Status Status;
	double* v =  (double*)malloc(N_part*sizeof(double));
	if(*rank == 0)
	{
		int dest;
		long long j, move = 1;
		for(dest=1; dest<=size; ++dest) // The last iteration is not sent to rank==size, but stays in rank==0.
		{
			#pragma omp parallel for schedule(static)
			for(j=0;j<N_part;++j)
				v[j] = 1.0/((double)(j+move)*(j+move));
			if(dest < size)
				MPI_Send(v,N_part,MPI_DOUBLE,dest,100,MPI_COMM_WORLD);
			move += N_part;
		}
	}else{
		MPI_Recv(v,N_part,MPI_DOUBLE,0,100,MPI_COMM_WORLD,&Status);
	}
	return v;
}

