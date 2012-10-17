
#include <R.h>

#include <stdio.h>
#include <stdlib.h>

void shooting(double* theta, double* y, double* x, int *m, int *n)
{
	int M, N;
	int i,j;
	double **X;

	M = *m;
	N = *n;

	// allocate and copy data matrix
	X = (double **) malloc(M * sizeof(double*));
	for(i = 0; i < M; i++) { 
		X[i] = (int *) malloc(N * sizeof(double)); 
		for( j=0; j<N; j++ ) {
			X[i][j] = x[ i + M * j ];
		}
	}
	
	
	for(i = 0; i < M; i++) {
		y[i] = 0;
		for( j=0; j< N; j++ ) {
			y[i] += X[i][j]*theta[j];
		}
	}

	// clean up	
	for(i = 0; i < N; i++){ free(X[i]); }
	free(X);
}

