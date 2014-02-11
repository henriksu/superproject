#include <stdio.h> // printf.
#include <math.h> // pow()
#include <stdlib.h> // malloc.

#include "sum.h" // All our functions + common.h ( MPI, omp, WallTime() ).

int main(int argc,char** argv){
	int i; // Generic loop variable.
	int N; // Numer of elements to sum.
	int method; // Indicates what summing medthod to use. 1=task 1, 2=task 2, 3=task 3, 4=task 4, 5=slowSum.
	int rank = 0; // MPI rank.
	double time; // Used to store both start time (seconds since 1970?) and time deltas in seconds.
	double S = (M_PI*M_PI)/6; // The limit of the series.
	double Sn; // Storing the partial sum.

	MPI_Init(&argc,&argv);

	// Parse commandline input.
	if(argc >= 2)
	{
		if (atoi(argv[1])>0 && atoi(argv[1])<6)
			method = atoi(argv[1]);
		else // Default
			method = 5;
		N = atoi(argv[2]);
	}
	else if (argc >=1)
	{
		method = 5;
		N = atoi(argv[1]);
	}		

	// Do calculations.
	switch (method)
	{
		case 1:
			printf("Running the non-parallelized programm (Task1).\n");
			time  = WallTime();
			Sn = sum(N);
			time = WallTime() - time;
			break;
		case 2:
			printf("Running the openMP-parallelized programm (Task2).\n");
			time= WallTime();
			Sn = sumShared(N);
			time = WallTime() - time;
			break;
		case 3:
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			if(rank==0)			
				printf("Running the MPI-parallelized programm (Task3).\n");
			time = WallTime();
			Sn = sumDist(N,&rank);
			time = WallTime() - time;
			break;
		case 4:
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			if(rank==0)
				printf("Running the openMP- and MPI-parallelized programm (Task4).\n");
			time = WallTime();
			Sn = sumHybrid(N,&rank);
			time = WallTime() - time;
			break;
		case 5:
			printf("Running the non-parallelized programm with better summation order(Task1).\n");
			time= WallTime();
			Sn = sumSlow(N);
			time = WallTime() - time;
			break;
	}

	// Output results.
	if(rank==0)
	{
		printf("n \terror \t\ttime\n");
		printf("%d \t%e \t%e\n", N, S-Sn, time);
	}
	
	// Clean up.
	MPI_Finalize();
	return 0;
}

