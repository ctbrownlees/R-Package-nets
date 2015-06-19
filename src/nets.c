
#include <R.h>
#include <math.h>

#define MAX(a,b)          (a>b?a:b)
#define MIN(a,b)          (a<b?a:b)
#define RHOIDX(i,j)       ( MAX(i,j)*(MAX(i,j)-1)/2 + MIN(i,j) )
#define ALPIDX(i,j,p,N,P) ( p*N*N + i*N + j )

double soft_thresholding(double c_yx,double c_xx,double lambda){
    double theta;
    if( c_yx > lambda/2 ){
        theta = (lambda/2 - c_yx)/c_xx;
    }
    else if( c_yx < -lambda/2 ){
        theta = (-lambda/2 - c_yx)/c_xx;
    }
    else{ // abs(c_yx) <= lambda/2
        theta = 0;
    }
    return theta;
}

//
void nets(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *kk, int *v)
{
	// variables 
	int T, N, P;
	double **y;
	double *y_aux, *x_aux;
	double *alpha_old;
	double *rho_old;
	double delta, toll;
	double lambda;
	int iter, maxiter;
	int verbose;
	int t, i, j, k, p;
	double rss, pen;
	double c_yx, c_xx;
  
	// init
	maxiter  = 100;
	toll     = 1e-4;
	verbose  = *v;
	T = *_T;
	N = *_N;
 	P = *_P;
	lambda = *_lambda;

	// allocate memory
	y = Calloc(T,double*);
	for( t=0; t<T; t++) { 
		y[t] = Calloc(N,double); 
		for( i=0; i<N; i++ ) {
			y[t][i] = _y[ T * i + t ];
		}
	}

	alpha_old = Calloc(N*N*P,double);
	rho_old   = Calloc(N*(N-1)/2,double);
	y_aux     = Calloc(T*N,double);
	x_aux     = Calloc(T*N,double);

	// Check A
	//for(p=0;p<P;++p){
	//	for(i=0;i<N;++i){
	//		for(j=0;j<N;++j){
	//			Rprintf( "%f, " , alpha[ ALPIDX(i,j,p,N,P) ]  );
	//		}
	//		Rprintf( "\n" );	
	//	}
	//	Rprintf( "\n" );
	//}
	//for( i=0;i<N*N*P; ++i) Rprintf( "%f, " , alpha[ i ]  );
	//Rprintf( "\n" );

	// Check R
	//for(i=0;i<N;++i){
	//	for(j=0;j<i;++j){
	//		Rprintf( "%f, " , rho[ RHOIDX(i,j) ]  );
	//	}
	//	Rprintf( "\n" );	
	//}
	//for( i=0;i<N*(N-1)/2; ++i) Rprintf( "%f, " , rho[ i ]  );

	// create residual
	memset(y_aux,0,sizeof(double)*N*T);
	for( i=0; i<N; ++i ){
		for( t=P; t<T; ++t ){

			y_aux[i*T+t] = y[t][i];

			// ALPHA
			for( j=0; j<N; ++j ){
				for( p=0; p<P; ++p ){
				  y_aux[i*T+t] -= alpha[ ALPIDX(i,j,p,N,P) ] * y[t-p-1][j];
				}
			}

			// RHO
			for( j=0; j<N; ++j){
				if( j!=i ){
					y_aux[i*T+t] -= rho[ RHOIDX(i,j) ] * sqrt(kk[i]/kk[j]) * y[t][j];
					for( p=0; p<P; ++p ){
						for( k=0; k<N; ++k ){
							y_aux[i*T+t] += rho[ RHOIDX(i,j) ] * sqrt(kk[i]/kk[j]) * alpha[ ALPIDX(j,k,p,N,P) ] * y[t-p-1][k];
						}
					}
				}
			}
		}
	}

	// main loop
	delta = 0.0;
	for( iter=1; iter<=maxiter; ++iter ){

		memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
		memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

		if( verbose ){

			rss = 0;
			pen = 0;
			for( i=0; i<N*T; ++i )       rss += (y_aux[i])*(y_aux[i]);
			for( i=0; i<N*N*P; ++i )     pen += lambda * alpha_weights[i] * fabs( alpha[i] );
			for( i=0; i<N*(N-1)/2; ++i ) pen += lambda * rho_weights[i]   * fabs( rho[i] )  ;

			Rprintf("Iter: %3.3d Rss %3.3f Pen %3.3f Obj %3.3f Delta: %f\n",iter,rss/(T*N),pen/(T*N),(rss+pen)/(T*N),delta);
		}
 
		// ALPHA 
		for( i=0; i<N; ++i){
			for( j=0; j<N; ++j){
				for( p=0; p<P; ++p ){

					// compute: y_aux, x_aux, c_yx, c_xx 
					c_yx = 0;
					c_xx = 0;
					for( k=0; k<N; ++k ){
						for( t=P; t<T; ++t ){
							if( k==i ){
								x_aux[ k*T+t ] = y[t-p-1][j];
							}
							else{
								x_aux[ k*T+t ] = -rho[ RHOIDX(i,k) ] * sqrt(kk[k]/kk[i]) * y[t-p-1][j];
							}
	
							y_aux[ k*T+t ] += alpha[ ALPIDX(i,j,p,N,P) ]*x_aux[ k*T+t ];
	
							c_yx -= y_aux[ k*T+t ] * x_aux[ k*T+t ];
							c_xx += x_aux[ k*T+t ] * x_aux[ k*T+t ];
						}
					}

					// update alpha
					alpha[ ALPIDX(i,j,p,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,p,N,P)]);
					//Rprintf("%d,%d,%d > c_yx %f c_xx %f lambda %f > alpha %f\n",i,j,p,c_yx,c_xx,lambda,alpha[ ALPIDX(i,j,p,N,P) ]);

					// update y_aux (if new coeff is != 0) 
					if( alpha[ ALPIDX(i,j,p,N,P) ] != 0.0 )
					{
						for( k=0; k<N; ++k ){
							for( t=P; t<T; ++t ){          
								y_aux[ k*T+t ] -= alpha[ ALPIDX(i,j,p,N,P) ]*x_aux[ k*T+t ];
							}
						}
					}

				}
			}
		}

		// RHO
		for( i=0; i<N; ++i){
			for( j=0; j<i; ++j ){

				// compute: y_aux, x_aux, c_yx, c_xx 
				c_yx = 0;
				c_xx = 0;
	
				memset(x_aux,0,sizeof(double)*N*T);				
				for( t=P; t<T; ++t ){
					
					x_aux[ i*T+t ] = y[t][j];
					x_aux[ j*T+t ] = y[t][i];

					// resid
					for( k=0; k<N; ++k ){
						for( p=0; p<P; ++p ){
							x_aux[ i*T+t ] -= alpha[ ALPIDX(j,k,p,N,P) ] * y[t-p-1][k];
							x_aux[ j*T+t ] -= alpha[ ALPIDX(i,k,p,N,P) ] * y[t-p-1][k];
						}
					}
					x_aux[ i*T+t ] = sqrt(kk[j]/kk[i]) * x_aux[ i*T+t ];
					x_aux[ j*T+t ] = sqrt(kk[i]/kk[j]) * x_aux[ j*T+t ];

					y_aux[ i*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ i*T+t ];
					y_aux[ j*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ j*T+t ];
	
					c_yx -= y_aux[ i*T+t ] * x_aux[ i*T+t ] + y_aux[ j*T+t ] * x_aux[ j*T+t ];
					c_xx += x_aux[ i*T+t ] * x_aux[ i*T+t ] + x_aux[ j*T+t ] * x_aux[ j*T+t ];
				}

				// update alpha
				rho[ RHOIDX(i,j) ] = soft_thresholding(c_yx,c_xx,lambda*rho_weights[RHOIDX(i,j)] ); 
				//Rprintf("%d,%d > c_yx %f c_xx %f lambda %f > rho %f\n",i,j,c_yx,c_xx,lambda,rho[RHOIDX(i,j)]);

				// update y_aux (if new coeff is != 0) 
				if( rho[ RHOIDX(i,j) ] != 0.0 )
				{
					for( t=P; t<T; ++t ){          
						y_aux[ i*T+t ] -= rho[ RHOIDX(i,j) ] * x_aux[ i*T+t ];
						y_aux[ j*T+t ] -= rho[ RHOIDX(i,j) ] * x_aux[ j*T+t ];
					}
				}

			}
		}
    
		// check for convergence
		delta = 0;
		for( i=0; i<N*N*P ; ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
		for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]); }
		if( delta<toll ) break;
	}
  
	// clean up	
	for(t = 0; t < T; t++){ Free(y[t]); }
	Free(y);
	Free(alpha_old);
	Free(rho_old);
	Free(y_aux);
	Free(x_aux);
	
}


/////////// OLD STUFF
char prog_indicator[] = {'-','\\','|','/','-','\\','|','/'};

// utility
double sgn(double x) {
    return (0 < x) - (x < 0);
}

// shooting algorithm
void shooting(double *theta, double *y, double *x, double *l, double *w, double *theta_init, int *m, int *n, int *v, int *ti, int *max_iter)
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
	maxiter = *max_iter;
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
			//Rprintf("X %d %d = %f. \n",i,j,X[i][j]);
		}
	}

	X_sumsqr = Calloc(N,double);
	for( j=0; j<N; j++) {
		X_sumsqr[j] = 0.0;
		for( i=0; i<M; i++ ) X_sumsqr[j] += X[i][j]*X[i][j];
		X_sumsqr[j] *= 2.0;
		//Rprintf("Xsquared %d = %f. \n",j,X_sumsqr[j]);
	}
	
	theta_cur = Calloc(N,double);

	// init
	for( j=0; j<N; ++j )
	{
		sum = 0.0;
		for( i=0; i<M; ++i ) sum += y[i]*X[i][j];
		sum *= -2.0;

		theta[j] = 0.0;

		if( sum - lambda*w[j] > 0 ) theta[j] = ( lambda*w[j]-sum)/X_sumsqr[j];
		if( sum + lambda*w[j] < 0 ) theta[j] = (-lambda*w[j]-sum)/X_sumsqr[j];
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
				for( k=0; k<N; ++k) {
					eps -= (k==j?0.0:theta[k])*X[i][k];
					//Rprintf("step %d %d %d eps :%f.\n",j,i,k,eps);
				}
				sum += eps*X[i][j];
				//Rprintf("step %d %d sum :%f.\n",j,i,sum);
			}
			sum *= -2.0;
			//Rprintf("step %d %d sum :%f,%f.\n",iter,j,abs(sum),sum);
			
			theta[j] = 0.0;
			if( sum - lambda*w[j] > 0 ) theta[j] = ( lambda*w[j]-sum)/X_sumsqr[j];
			if( sum + lambda*w[j] < 0 ) theta[j] = (-lambda*w[j]-sum)/X_sumsqr[j];
			//Rprintf("theta was %f now it is %f. \n",theta_cur[j],theta[j]);
		}

		// check for convergence
		delta = 0.0;
		for( j=0; j<N; ++j ) delta += (theta_cur[j]-theta[j])*(theta_cur[j]-theta[j]); 
		delta = sqrt( delta );
		if( delta < toll ) break;
	}
	Rprintf("\bdone!");

	if( maxiter==iter ) {
		Rprintf("Convergence was not achieved after %d iterations.\n",maxiter);
	}

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
	maxiter = 200;
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
					for(k=0;k<N;++k) eps -= (k==j?0.0:theta[k])*X[i][k];
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
				for(k=0;k<N;++k) eps -= (k==j?0.0:theta[k])*X[i][k];
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
	double *theta_act;
	double toll;
	double lambda;
	double delta;
	int iter_in, maxiter_in, iter_out, maxiter_out;
	int verbose;
	int i, j, k, h;
	double m1, m2;
	double cov, var;
	double eps, eps_1, eps_2;

	int n_active;

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
	theta_act = Calloc(N*(N-1)/2,double);

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

		n_active = 0;
		for( i=1; i<N; ++i ) {
			for( j=0; j<i-1; ++j ) {
				theta_cur[ i*(i-1)/2+j ] = theta[ i*(i-1)/2+j ];
				theta_act[i*(i-1)/2+j ] = (double)(theta[i*(i-1)/2+j ]!=0);
				n_active    += (int) (theta[i*(i-1)/2+j ]!=0);
			}
		}		

		if(verbose) Rprintf(" > outer iter: %d %d active parameters out of %d\n",iter_out+1,n_active,N);

		// inner theta loop
		for( iter_in=0 ; iter_in<maxiter_in ; ++iter_in ) {

			// set current param
			for( i=1; i<N; ++i ) for( j=0; j<i-1; ++j ) theta_cur[ i*(i-1)/2+j ] = theta[ i*(i-1)/2+j ];

			// Active Set UPDATE
			for( i=1; i<N; ++i ){
				for( j=0; j<i-1; ++j ){

					if( theta_act[i*(i-1)/2+j]==0 ) continue;

					// nasty bit!
					cov = 0.0;
					var = 0.0;	
					for( k=0; k<M; ++k ){
						eps_1 = Y[k][i];
						eps_2 = Y[k][j];
						for( h=0; h<N; ++h ){
							if( h==j || h==i ) continue;
							eps_1 -= theta[ (h<i)?(i*(i-1)/2+h):(h*(h-1)/2+i) ] * sqrt(ivar[h]/ivar[i]) * (Y[k][h]);
							eps_2 -= theta[ (h<j)?(j*(j-1)/2+h):(h*(h-1)/2+j) ] * sqrt(ivar[h]/ivar[j]) * (Y[k][h]);
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

		// set current param
		for( i=1; i<N; ++i ) for( j=0; j<i-1; ++j ) theta_cur[ i*(i-1)/2+j ] = theta[ i*(i-1)/2+j ];

		// Non Active Set UPDATE
		for( i=1; i<N; ++i ){
			for( j=0; j<i-1; ++j ){

				if( theta_act[i*(i-1)/2+j]!=0 ) continue;

				// nasty bit!
				cov = 0.0;
				var = 0.0;	
				for( k=0; k<M; ++k ){
					eps_1 = Y[k][i];
					eps_2 = Y[k][j];
					for( h=0; h<N; ++h ){
						if( h==j || h==i ) continue;
						eps_1 -= theta[ (h<i)?(i*(i-1)/2+h):(h*(h-1)/2+i) ] * sqrt(ivar[h]/ivar[i]) * (Y[k][h]);
						eps_2 -= theta[ (h<j)?(j*(j-1)/2+h):(h*(h-1)/2+j) ] * sqrt(ivar[h]/ivar[j]) * (Y[k][h]);
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

		// update variances
		for( i=0; i<N; ++i ){
			var = 0.0;	
			for( k=0; k<M; ++k ){
				eps_1 = Y[k][i];
				for( j=0; j<N; ++j ){
					if( i==j ) continue;
					eps_1 -= theta_cur[ (j<i)?(i*(i-1)/2+j):(j*(j-1)/2+i) ] * sqrt(ivar[j]/ivar[i]) * ( Y[k][j] );
				}
				var += eps_1*eps_1;
			}
			var /= ((double) M);
			ivar[i] = 1.0 / var;
		}

		// check for convergence
		delta = 0.0;
		for( i=1; i<N; ++i ) for( j=0; j<i-1; ++j ) delta += (theta_cur[i*(i-1)/2+j]-theta[i*(i-1)/2+j])*(theta_cur[i*(i-1)/2+j]-theta[i*(i-1)/2+j]); 
		delta = sqrt( delta );	
		if(verbose) Rprintf("iter: %3d toll: %f\n",iter_in+1,delta);
		if( delta < toll ) break;

	}

	// clean up	
	for(i = 0; i < M; i++){ Free(Y[i]); }
	Free(Y);
	Free(theta_cur);
}

