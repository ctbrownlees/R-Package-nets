
#include <R.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void shooting(double *theta, double *y, double *x, double *l, double *w, int *m, int *n)
{
	int M, N;
	int i,j;
	double **X;
	double *X_sumsqr;
	double *theta_cur;
	double lambda;
	int iter, maxiter;
	double toll, delta;
	double sum, eps;

	// 
	maxiter = 100;
	toll = 1e-6;

	// housekeeping
	M = *m;
	N = *n;
	lambda = *l;

	// allocate stuff
	X = (double **) malloc(M * sizeof(double*));
	for( i=0; i<M; i++) { 
		X[i] = (double *) malloc(N * sizeof(double)); 
		for( j=0; j<N; j++ ) {
			X[i][j] = x[ M * j + i ];
		}
	}

	X_sumsqr = (double *) malloc(N * sizeof(double));
	for( j=0; j<N; j++) {
		X_sumsqr[j] = 0.0;
		for( i=0; i<M; i++ ) X_sumsqr[j] += X[i][j]*X[i][j];
		X_sumsqr[j] *= 2;
	}
	
	theta_cur = (double *) malloc(N * sizeof(double));

	// init
	for( j=0; j<N; ++j )
	{
		sum = 0.0;
		for( i=0; i<M; ++i ) sum += y[i]*X[i][j];
		sum *= -2;

		if( abs(sum)>lambda*w[j] )
		{
			if( sum - lambda*w[j] > 0 ) theta[j] = ( lambda*w[j]-sum)/X_sumsqr[j];
			if( sum + lambda*w[j] < 0 ) theta[j] = (-lambda*w[j]-sum)/X_sumsqr[j];
		}
		else theta[j] = 0.0;
	}

	// iterate
	for( iter=0 ; iter<maxiter ; ++iter )
	{
		memcpy( theta_cur , theta , N*sizeof(double) );

		sum = 0.0;
		for( i=0; i<M; ++i )
		{
			eps = y[i];
			for(j=0;j<N;++j) eps -= (i==j?theta_cur[i]:0.0)*X[i][j];
			sum += eps*X[i][j];
		}
		sum *= -2;

		if( abs(sum)>lambda*w[j] )
		{
			if( sum - lambda*w[i] > 0 ) theta[j] = ( lambda*w[i]-sum)/X_sumsqr[j];
			if( sum + lambda*w[i] < 0 ) theta[j] = (-lambda*w[i]-sum)/X_sumsqr[j];
		}
		else theta[j] = 0.0;

		// check for convergence
		delta = 0.0;
		for( j=0; j<N; ++j ) delta += (theta_cur[j]-theta[j])*(theta_cur[j]-theta[j]); 
		delta = sqrt( delta );
	
		if( delta < toll ) break;

		printf("iter: %3d toll: %f\n",iter+1,delta);
	}

	// clean up	
	for(i = 0; i < N; i++){ free(X[i]); }
	free(X);
	free(X_sumsqr);
	free(theta_cur);
}

