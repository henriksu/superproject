#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc,char** argv){

	int N = pow(2,14);
	int indexes[11];
	int i;
	double S = (M_PI*M_PI)/6;
	double* v = (double*)malloc(N*sizeof(double));
	
	indexes[0] = 8;
	for(i=1;i<11;++i){
		indexes[i] = 2*indexes[i-1];
	}
	for(i=0;i<N;++i){
		v[i] = 1.0/(i*i);
	}
	double* Sn = (double*)calloc(11,sizeof(double));
	int j=10;
	for(i=N-1;i>=0;--i){
		Sn[j] += v[i];
		if(i==indexes[j-1]){
			--j;
		}
	}
	for(i=1;i<11;++i){
		Sn[i] += Sn[i-1];
	}
	

	for(i=0;i<11;++i){
		printf("%f\n",S-Sn[i]);
	}
	return 0;
}
