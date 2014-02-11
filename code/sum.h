#include "common.h" // MPI, omp, WallTime().

// Implemented in sum.c.
// All these functions return the partial sum of the first N elements of the seqence 1/(i*i).
double sum(const long long N);
double sumSlow(const long long N);
double sumShared(const long long N);
double sumDist(const long long N, int* rank);
double sumHybrid(const long long N, int* rank);

// Implemented in initVec.c.
// All these functions return a pointer to an array of length N or N_part.
double* initVec_seq(const long long N);
double* initVec_shared(const long long N);
double* initVec_dist(const long long N_part, const int* rank, const int size);
double* initVec_hybrid(const long long N_part, const int* rank, const int size);

// Implemented in summing.c.
// All these functions return the sum of the elements in the array.
double summingSeq(double* v, const long long N);
double summingSlow(double* v, const long long N);
double summingShared(double* v, const long long N);
