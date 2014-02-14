#include <stdio.h> // printf.
#include <math.h> // pow()
#include <stdlib.h> // malloc.

#include "sum.h" // All our functions + common.h ( MPI, omp, WallTime() ).

int main(int argc,char** argv){
	int iterations = 30;
	long long N[iterations]; // Numer of elements to sum.
	int i;
	int rank = 0; // MPI rank.
	double S = (M_PI*M_PI)/6; // The limit of the series.
	
	double* SnOrdi = (double*)malloc(iterations*sizeof(double));
	double* SnForw = (double*)malloc(iterations*sizeof(double));

	MPI_Init(&argc,&argv);

	N[0] = 8;
	for(i=1; i<iterations; ++i)
		N[i] = 2*N[i-1];

	// Do calculations.
	for(i=0; i<iterations; ++i)
	{
		SnOrdi[i] = sum(N[i]);
		SnForw[i] = sumF(N[i]);
	}
	// Output results.

	if(rank==0)
	{
		printf("n \terror\n");
		for(i=0; i<iterations; ++i)
			printf("%lld \t%.17e \t%.17e \n", N[i], S-SnOrdi[i], S-SnForw[i]);
	}
	
	// Clean up.
	free(SnOrdi); free(SnForw);
	MPI_Finalize();
	return 0;
}

