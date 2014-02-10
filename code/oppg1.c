#include <stdio.h> // printf.
#include <math.h> // pow()
#include <stdlib.h> // malloc.

#include "sum.h" // All our functions + common.h ( MPI, omp, WallTime() ).
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





