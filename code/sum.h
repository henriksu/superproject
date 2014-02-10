#include "common.h" // MPI, omp, WallTime().

// Implemented in sum.c.
// All these functions return the partial sum of the first N elements of the seqence 1/(i*i).
double sum(const int N);
double sumSlow(const int N);
double sumShared(const int N);
double sumDist(const int N, int* rank);
double sumHybrid(const int N, int* rank);

// Implemented in initVec.c.
// All these functions return a pointer to an array of length N or N_part.
double* initVec_seq(const int N);
double* initVec_shared(const int N);
double* initVec_dist(const int N_part, const int* rank, const int size);
double* initVec_hybrid(const int N_part, const int* rank, const int size);

// Implemented in summing.c.
// All these functions return the sum of the elements in the array.
double summingSeq(double* v, const int N);
double summingSlow(double* v, const int N);
double summingShared(double* v, const int N);
