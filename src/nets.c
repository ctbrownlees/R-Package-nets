
#include <R.h>

#include <math.h>

char prog_indicator[] = {'-','\\','|','/','-','\\','|','/'};

// utility
double sgn(double x) {
    return (0 < x) - (x < 0);
}

// shooting algorithm
void shooting(double *theta, double *y, double *x, double *l, double *w, double *theta_init, int *m, int *n, int *v, int *ti)
{
	int M, N;
	int i,j,k;
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
	X = Calloc(M,double*);
	for( i=0; i<M; i++) { 
		X[i] = Calloc(N,double); 
		for( j=0; j<N; j++ ) {
			X[i][j] = x[ M * j + i ];
		}
	}

	X_sumsqr = Calloc(N,double);
	for( j=0; j<N; j++) {
		X_sumsqr[j] = 0.0;
		for( i=0; i<M; i++ ) X_sumsqr[j] += X[i][j]*X[i][j];
		X_sumsqr[j] *= 2;
	}
	
	theta_cur = Calloc(N,double);

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
	Rprintf("\r LASSO Optimization (shooting algorithm):  ");
	for( iter=0 ; iter<maxiter ; ++iter ) {
		Rprintf("\b%c",prog_indicator[ iter % 7 ]);

		for( j=0; j<N; ++j ) theta_cur[j] = theta[j];

		for( j=0; j<N; ++j ) {
			sum = 0.0;
			for( i=0; i<M; ++i ) {
				eps = y[i];
				for(k=0;k<N;++k) eps -= (k==j?0.0:theta_cur[k])*X[i][k];
				sum += eps*X[i][j];
			}
			sum *= -2;

			if( abs(sum)>lambda*w[j] ) {
				if( sum - lambda*w[j] > 0 ) theta[j] = ( lambda*w[j]-sum)/X_sumsqr[j];
				if( sum + lambda*w[j] < 0 ) theta[j] = (-lambda*w[j]-sum)/X_sumsqr[j];
			}
			else theta[j] = 0.0;
		}

		// check for convergence
		delta = 0.0;
		for( j=0; j<N; ++j ) delta += (theta_cur[j]-theta[j])*(theta_cur[j]-theta[j]); 
		delta = sqrt( delta );
		if( delta < toll ) break;
	}
	Rprintf("\bdone!");

	// clean up	
	for(i = 0; i < M; i++){ Free(X[i]); }
	Free(X);
	Free(X_sumsqr);
	Free(theta_cur);
}

void activeshooting(double *theta, double *y, double *x, double *l, double *w, double *theta_init, int *m, int *n, int *v, int *ti)
{
	int M, N;
	int i,j,k;
	double **X;
	double *X_sumsqr;
	double *theta_cur;
	double *theta_act;
	double lambda;
	int verbose;
	int out_iter, in_iter, maxiter;
	double toll, delta;
	double sum, eps;

	int n_active;

	// 
	maxiter = 100;
	toll = 1e-6;
	verbose = *v;

	// housekeeping
	M = *m;
	N = *n;
	lambda = *l;

	if( verbose ) Rprintf("\nActive Shooting Algorithm!\n");

	// allocate stuff
	X = Calloc(M,double*);
	for( i=0; i<M; i++) { 
		X[i] = Calloc(N,double); 
		for( j=0; j<N; j++ ) {
			X[i][j] = x[ M * j + i ];
		}
	}

	X_sumsqr = Calloc(N,double);
	for( j=0; j<N; j++) {
		X_sumsqr[j] = 0.0;
		for( i=0; i<M; i++ ) X_sumsqr[j] += X[i][j]*X[i][j];
		X_sumsqr[j] *= 2;
	}
	
	theta_cur = Calloc(N,double);
	theta_act = Calloc(N,double);

	// init
	n_active = 0;
	if( verbose ) Rprintf(" > initialising parameters");
	for( j=0; j<N; ++j )
	{
		if( verbose ) if( (j % 10)==0 ) Rprintf(".");

		sum = 0.0;
		for( i=0; i<M; ++i ) sum += y[i]*X[i][j];
		sum *= -2;

		if( abs(sum)>lambda*w[j] ) {
			if( sum - lambda*w[j] > 0 ) theta[j] = ( lambda*w[j]-sum)/X_sumsqr[j];
			if( sum + lambda*w[j] < 0 ) theta[j] = (-lambda*w[j]-sum)/X_sumsqr[j];
		}
		else theta[j] = 0.0;
	}
	if( verbose ) Rprintf("\n");

	// iterate
	for( out_iter=0 ; out_iter<maxiter ; ++out_iter ) {
	
		// 
		n_active = 0;
		for( j=0; j<N; ++j ){ 
			theta_cur[j] = theta[j];
			theta_act[j] = (double)(theta[j]!=0);
			n_active    += (int) (theta[j]!=0);
		}

		if(verbose) Rprintf(" > outer iter: %d %d active parameters out of %d\n",out_iter+1,n_active,N);

		// Active Set UPDATE
		for( in_iter=0 ; in_iter<maxiter ; ++in_iter ) {

			for( j=0; j<N; ++j ) theta_cur[j] = theta[j];

			for( j=0; j<N; ++j ) {
				if( theta_act[j]==0 ) continue;

				sum = 0.0;
				for( i=0; i<M; ++i ) {
					eps = y[i];
					for(k=0;k<N;++k) eps -= (k==j?0.0:theta_cur[k])*X[i][k];
					sum += eps*X[i][j];
				}
				sum *= -2;

				if( abs(sum)>2*lambda*w[j] ) {
					theta[j] = sgn(-sum) * ( abs(sum)/X_sumsqr[j] - (2.0*lambda*w[j])/X_sumsqr[j] );
				}
				//if( abs(sum)>lambda*w[j] ) {
					//if( sum - lambda*w[j] > 0 ) theta[j] = ( lambda*w[j]-sum)/X_sumsqr[j];
					//if( sum + lambda*w[j] < 0 ) theta[j] = (-lambda*w[j]-sum)/X_sumsqr[j];
				//}	
				else theta[j] = 0.0;
			}
	
			delta = 0.0;
			for( j=0; j<N; ++j ) delta += (theta_cur[j]-theta[j])*(theta_cur[j]-theta[j]); 
			delta = sqrt( delta );	
			if( delta < toll ) break;
		}

		// Non Active Set UPDATE
		for( j=0; j<N; ++j ) theta_cur[j] = theta[j];

		for( j=0; j<N; ++j ) {
			if( theta_act[j]!=0 ) continue;

			sum = 0.0;
			for( i=0; i<M; ++i )
			{
				eps = y[i];
				for(k=0;k<N;++k) eps -= (k==j?0.0:theta_cur[k])*X[i][k];
				sum += eps*X[i][j];
			}
			sum *= -2;

			//if( abs(sum)>lambda*w[j] )
			//{
			//	if( sum - lambda*w[j] > 0 ) theta[j] = ( lambda*w[j]-sum)/X_sumsqr[j];
			//	if( sum + lambda*w[j] < 0 ) theta[j] = (-lambda*w[j]-sum)/X_sumsqr[j];
			//}
			if( abs(sum)>2*lambda*w[j] ) {
				theta[j] = sgn(-sum) * ( abs(sum)/X_sumsqr[j] - (2*lambda*w[j])/X_sumsqr[j] );
			}
			else theta[j] = 0.0;
		}

		// check for convergence
		delta = 0.0;
		for( j=0; j<N; ++j ) delta += (theta_cur[j]-theta[j])*(theta_cur[j]-theta[j]); 
		delta = sqrt( delta );	
		if(verbose) Rprintf(" > toll: %f\n",delta);

		if( delta < toll ) break;
	}

	// clean up	
	for(i = 0; i < M; i++){ Free(X[i]); }
	Free(X);
	Free(X_sumsqr);
	Free(theta_cur);
	Free(theta_act);
}

//
void space(double *theta, double *ivar, double *y, double *l, int *m, int *n, int *v)
{
	int M, N;
	double **Y;
	double *theta_cur;
	double toll;
	double lambda;
	double delta;
	int iter_in, maxiter_in, iter_out, maxiter_out;
	int verbose;
	int i, j, k, h;
	double m1, m2;
	double cov, var;
	double eps, eps_1, eps_2;

	// 
	maxiter_out = 5;
	maxiter_in  = 100;
	toll = 1e-6;

	// housekeeping
	verbose = *v;
	M = *m;
	N = *n;
	lambda = *l;

	if( verbose ) Rprintf("\nSPaCE Algorithm!\n");

	// allocate stuff
	Y = Calloc(M,double*);
	for( i=0; i<M; i++) { 
		Y[i] = Calloc(N,double); 
		for( j=0; j<N; j++ ) {
			Y[i][j] = y[ M * j + i ];
		}
	}
	theta_cur = Calloc(N*(N-1)/2,double);

	// init variances
	for( i=0; i<N; ++i ){
		m2=0;
		m1=0;
		for(k=0;k<M;++k){
			m2 += Y[k][i]*Y[k][i];
			m1 += Y[k][i];
		}
		m2 /= ((double) M);
		m1 /= ((double) M);
		ivar[i] = 1.0 / ( m2 - m1*m1 );
	}

	// init theta
	for( i=1; i<N; ++i )
	{
		for( j=0; j<i; ++j ){
			cov = 0.0;
			var = 0.0;
			for( k=0; k<M; ++k ){
				cov += (Y[k][i]) * sqrt(ivar[j]/ivar[i])*(Y[k][j]) + (Y[k][j]) * sqrt(ivar[i]/ivar[j])*(Y[k][i]);
				var += (ivar[j]/ivar[i])*(Y[k][j]*Y[k][j]) + (ivar[i]/ivar[j])*(Y[k][i]*Y[k][i]);
			}
			cov *= -2.0;
			var *= +2.0;

			if( abs(cov)>lambda ){
				if( cov - lambda > 0 ) theta[ i*(i-1)/2 + j ] = ( lambda-cov)/var;
				if( cov + lambda < 0 ) theta[ i*(i-1)/2 + j ] = (-lambda-cov)/var;
			}
			else theta[ i*(i-1)/2 +j ] = 0.0;
		}
	}

	// space outer loop
	for( iter_out=0 ; iter_out<maxiter_out ; ++iter_out ) {

		// inner theta loop
		for( iter_in=0 ; iter_in<maxiter_in ; ++iter_in ) {

			// set current param
			for( i=1; i<N; ++i ) for( j=0; j<i-1; ++j ) theta_cur[ i*(i-1)/2 +j ] = theta[ i*(i-1)/2+j ];

			// update
			for( i=1; i<N; ++i ){
				for( j=0; j<i-1; ++j ){

					// nasty bit!
					cov = 0.0;
					var = 0.0;	
					for( k=0; k<M; ++k ){
						eps_1 = Y[k][i];
						eps_2 = Y[k][j];
						for( h=0; h<N; ++h ){
							if( h==j || h==i ) continue;
							eps_1 -= theta_cur[ (h<i)?(i*(i-1)/2+h):(h*(h-1)/2+i) ] * sqrt(ivar[h]/ivar[i]) * (Y[k][h]);
							eps_2 -= theta_cur[ (h<j)?(j*(j-1)/2+h):(h*(h-1)/2+j) ] * sqrt(ivar[h]/ivar[j]) * (Y[k][h]);
						}
						cov += eps_1 * sqrt(ivar[j]/ivar[i]) * (Y[k][j]) + eps_2 * sqrt(ivar[i]/ivar[j]) * (Y[k][i]);
						var += (ivar[j]/ivar[i])*(Y[k][j]*Y[k][j]) + (ivar[i]/ivar[j])*(Y[k][i]*Y[k][i]);
					}
					cov *= -2.0;
					var *= +2.0;

					// shoot!
					if( abs(cov)>lambda ){
						if( cov - lambda > 0 ) theta[i*(i-1)/2+j] = ( lambda-cov)/var;
						if( cov + lambda < 0 ) theta[i*(i-1)/2+j] = (-lambda-cov)/var;
					}
					else theta[i*(i-1)/2+j] = 0.0;
				}
			}

			// check for convergence
			delta = 0.0;
			for( i=1; i<N; ++i ) for( j=0; j<i-1; ++j ) delta += (theta_cur[i*(i-1)/2+j]-theta[i*(i-1)/2+j])*(theta_cur[i*(i-1)/2+j]-theta[i*(i-1)/2+j]); 
			delta = sqrt( delta );	
			if(verbose) Rprintf("iter: %3d toll: %f\n",iter_in+1,delta);
			if( delta < toll ) break;
		}

		// update variances
		for( i=0; i<N; ++i ){
			var = 0.0;	
			for( k=0; k<M; ++k ){
				eps_1 = Y[k][i];
				for( j=0; j<N; ++j ){
					if( i==j ) continue;
					eps_1 -= theta_cur[ (j<i)?(i*(i-1)/2+j):(j*(j-1)/2+i) ] * sqrt(ivar[j]/ivar[i]) * (Y[k][j]);
				}
				var += eps_1*eps_1;
			}
			var /= ((double) M);
			ivar[i] = 1.0 / var;
		}

	}

	// clean up	
	for(i = 0; i < M; i++){ Free(Y[i]); }
	Free(Y);
	Free(theta_cur);
}

