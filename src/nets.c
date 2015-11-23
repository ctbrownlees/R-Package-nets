
#include <R.h>
#include <math.h>

#define MAX(a,b)          (a>b?a:b)
#define MIN(a,b)          (a<b?a:b)
#define RHOIDX(i,j)       ( MAX(i,j)*(MAX(i,j)-1)/2 + MIN(i,j) )
#define ALPIDX(i,j,p,N,P) ( p*N*N + i*N + j )

void nets_shooting(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v, int *m);
void nets_activeshooting(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v, int *m);

double soft_thresholding(double c_yx,double c_xx,double lambda);

void alpha_update(double *alpha, int i, int j, int k, double *y_aux, double *x_aux, double *rho, double *c, double **y, double lambda, double *alpha_weights, int T, int N, int P);
void rho_update(double *rho, int i, int j, double *y_aux, double *x_aux, double *alpha, double *c, double **y, double lambda, double *rho_weights, int T, int N, int P);

void nets_log(double *alpha, double *rho, double *y_aux, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int GN, int CN, int iter, double delta);
void nets_sanity_check(double **y, double *alpha, double *rho, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int GN, int CN);

// soft thresholding operator
double soft_thresholding(double c_yx,double c_xx,double lambda){
    double theta;
    if( -c_yx > lambda/2 ){
        theta = (lambda/2 + c_yx)/c_xx;
    }
    else if( -c_yx < -lambda/2 ){
        theta = (-lambda/2 + c_yx)/c_xx;
    }
    else{ // abs(c_yx) <= lambda/2
        theta = 0;
    }
    return theta;
}

// (simple) SHOOTING
void nets_shooting(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v, int *m) {

	// variables 
	int T, N, P;
	int granger_network, parcorr_network;
	double **y;
	double *y_aux, *x_aux;
	double *alpha_old;
	double *rho_old;
	double delta, toll;
	double lambda;
	int iter, maxiter;
	int verbose;
	int t, i, j, k, p;
	double c_yx, c_xx;
  
	// init
	maxiter  = *m;
	toll     = 1e-4;
	verbose  = *v;
	T = *_T;
	N = *_N;
 	P = *_P;
	lambda = *_lambda;
	granger_network = *GN;
	parcorr_network = *CN;

	// allocate memory
	y = (double **) R_alloc(T,sizeof(double*));
	for( t=0; t<T; t++) { 
		y[t] = (double *) R_alloc(N,sizeof(double));
		for( i=0; i<N; i++ ) {
			y[t][i] = _y[ T * i + t ];
		}
	}

	alpha_old = (double *) R_alloc(N*N*P,sizeof(double));
	rho_old   = (double *) R_alloc(N*(N-1)/2,sizeof(double));
	y_aux     = (double *) R_alloc(T*N,sizeof(double));
	x_aux     = (double *) R_alloc(T*N,sizeof(double));

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
					y_aux[i*T+t] -= rho[ RHOIDX(i,j) ] * sqrt(c[j]/c[i]) * y[t][j];
					for( p=0; p<P; ++p ){
						for( k=0; k<N; ++k ){
							y_aux[i*T+t] += rho[ RHOIDX(i,j) ] * sqrt(c[j]/c[i]) * alpha[ ALPIDX(j,k,p,N,P) ] * y[t-p-1][k];
						}
					}
				}
			}
		}
	}

	// main loop
	for( iter=1; iter<=maxiter; ++iter ){

		if( verbose ) nets_log(alpha,rho,y_aux,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,iter,delta);

		memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
		memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

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
								x_aux[ k*T+t ] = -rho[ RHOIDX(i,k) ] * sqrt(c[k]/c[i]) * y[t-p-1][j];
							}
	
							y_aux[ k*T+t ] += alpha[ ALPIDX(i,j,p,N,P) ]*x_aux[ k*T+t ];
	
							c_yx += y_aux[ k*T+t ] * x_aux[ k*T+t ];
							c_xx += x_aux[ k*T+t ] * x_aux[ k*T+t ];
						}
					}

					// update alpha
					alpha[ ALPIDX(i,j,p,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,p,N,P)]);

					//Rprintf("VERY EXACT: %d %d %d -> %f %f : beta_ls %f beta_lasso %f\n",1+i,1+j,1+p,c_yx,c_xx,c_yx/c_xx,alpha[ ALPIDX(i,j,p,N,P) ]);

					if( alpha[ ALPIDX(i,j,p,N,P) ] != 0.0 )
					{
						for( k=0; k<N; ++k ){
							for( t=P; t<T; ++t ){
								y_aux[ k*T+t ] -= alpha[ ALPIDX(i,j,p,N,P) ] * x_aux[ k*T+t ];
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
					x_aux[ i*T+t ] = sqrt(c[j]/c[i]) * x_aux[ i*T+t ];
					x_aux[ j*T+t ] = sqrt(c[i]/c[j]) * x_aux[ j*T+t ];

					y_aux[ i*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ i*T+t ];
					y_aux[ j*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ j*T+t ];
	
					c_yx += y_aux[ i*T+t ] * x_aux[ i*T+t ] + y_aux[ j*T+t ] * x_aux[ j*T+t ];
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

	if( verbose ) nets_log(alpha,rho,y_aux,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,0,delta);	
}

// ACTIVE SHOOTING
void nets_activeshooting(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v, int *m)
{
	// variables 
	int T, N, P;
	int granger_network, parcorr_network;
	double **y, **eps;
	double *alpha_old, *rho_old;
	double **alpha_set, **rho_set;
	int alpha_nact, rho_nact;
	double delta, toll;
	double lambda;
	int iter1, iter2, maxiter;
	int verbose;
	int idx, i, j, k, l;
	int t;

	double *y_aux;
	double *x_aux;
  
	// init
	maxiter  = *m;
	toll     = 1e-4;
	verbose  = *v;
	T = *_T;
	N = *_N;
 	P = *_P;
	lambda = *_lambda;
	granger_network = *GN;
	parcorr_network = *CN;

	if( granger_network==0 ){ P = 0; };

	// allocate memory
	y   = (double **) R_alloc(T,sizeof(double*));
	eps = (double **) R_alloc(T,sizeof(double*));
	for( t=0; t<T; t++) { 
		y[t]   = (double *) R_alloc(N,sizeof(double)); 
		eps[t] = (double *)R_alloc(N,sizeof(double)); 
		for( i=0; i<N; i++ ) {
			y[t][i]   = _y[ T * i + t ];
			eps[t][i] = 0;
		}
	}

	y_aux = (double *) R_alloc(T*N,sizeof(double));
	memset(y_aux,0,sizeof(double)*N*T);				
	x_aux = (double *) R_alloc(T*N,sizeof(double));
	
	if( granger_network ){
		alpha_old = (double *) R_alloc(N*N*P,sizeof(double));

		alpha_set = (double **) R_alloc(N*N*P,sizeof(double *));
		for( idx=0; idx<N*N*P; ++idx){
			alpha_set[idx] = (double *) R_alloc(3,sizeof(double)); 
		}
	}
	
	if( parcorr_network ){
		rho_old   = (double *) R_alloc(N*(N-1)/2,sizeof(double));

		rho_set = (double **) R_alloc(N*(N-1)/2,sizeof(double *));
		for( idx=0; idx<N*(N-1)/2; ++idx){
			rho_set[idx] = (double *) R_alloc(2,sizeof(double)); 
		}
	}

	// check! 
	//nets_sanity_check(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network);

	// ALPHA INIT
	/*
	if( granger_network ){
		for( i=0; i<N; ++i) for( j=0; j<N; ++j) for( k=0; k<P; ++k ) 1+1;// alpha_update_init(alpha,i,j,k,y,lambda,alpha_weights,T,N,P);
	}

	// RHO INIT
	if( parcorr_network ){
		for( i=0; i<N; ++i) for( j=0; j<i; ++j) 1+1;//rho_update_init(rho,i,j,y_aux,x_aux,alpha,c,y,lambda,rho_weights,T,N,P);
	}
	*/

	// y_aux INIT
	for( i=0; i<N; ++i ){
		for( t=P; t<T; ++t ){

			y_aux[i*T+t] = y[t][i];

			// ALPHA
			if( granger_network ){
				for( j=0; j<N; ++j ){
					for( k=0; k<P; ++k ){
					  y_aux[i*T+t] -= alpha[ ALPIDX(i,j,k,N,P) ] * y[t-k-1][j];
					}
				}
			}

			// RHO
			if( parcorr_network ){
				for( k=0; k<P; ++k ){
					for( j=0; j<N; ++j){
						for( l=0; l<N; ++l ){
							if( i!=l ) y_aux[i*T+t] += rho[ RHOIDX(i,l) ] * sqrt(c[l]/c[i]) * alpha[ ALPIDX(l,j,k,N,P) ] * y[t-k-1][j];
						}
					}
				}
				for( l=0; l<N; ++l ){
					if( i!=l ) y_aux[i*T+t] -= rho[ RHOIDX(i,l) ] * sqrt(c[l]/c[i]) * y[t][l];
				}
			}
		}
	}

	// Main Loop
	for( iter1=1; iter1<=maxiter; ++iter1 ){
			
		if( verbose ) nets_log(alpha,rho,y_aux,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,iter1,delta);

		// ACTIVE Set 
		if( granger_network ){
			alpha_nact = 0;
			for( i=0; i<N; ++i){ 
				for( j=0; j<N; ++j){ 
					for( k=0; k<P; ++k ){
						if( alpha[ ALPIDX(i,j,k,N,P) ]!=0 ){
							alpha_set[ alpha_nact ][0] = i;
							alpha_set[ alpha_nact ][1] = j;
							alpha_set[ alpha_nact ][2] = k;
							++alpha_nact;
						}
					}
				}
			}
		}
		if( parcorr_network ){
			rho_nact = 0;
			for( i=0; i<N; ++i){ 
				for( j=0; j<i; ++j){ 
					if( rho[ RHOIDX(i,j) ]!=0 ){
						rho_set[ rho_nact ][0] = i;
						rho_set[ rho_nact ][1] = j;
						++rho_nact;
					}
				}
			}
		}

		// ACTIVE Update 
		for( iter2=1; iter2<=maxiter; ++iter2 ){

			if( granger_network ) memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
			if( parcorr_network ) memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

			// ALPHA Step
			if( granger_network ){
				for( idx=0; idx<alpha_nact; ++idx){
					i = alpha_set[idx][0];
					j = alpha_set[idx][1];
					k = alpha_set[idx][2];

					alpha_update(alpha,i,j,k,y_aux,x_aux,rho,c,y,lambda,alpha_weights,T,N,P);
				}
			}

			// RHO Step
			if( parcorr_network ){
				for( idx=0; idx<rho_nact; ++idx){
					i = rho_set[idx][0];
					j = rho_set[idx][1];

					rho_update(rho,i,j,y_aux,x_aux,alpha,c,y,lambda,rho_weights,T,N,P);
				}
			}

			// Convergence Check
			delta = 0;
			if( granger_network ) for( i=0; i<N*N*P ;    ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
			if( parcorr_network ) for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]);     }
			if( delta<toll ) break;
			if( iter1==1 ) break;
		}

		// GLOBAL UPDATE
		if( granger_network ) memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
		if( parcorr_network ) memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

		// ALPHA Step
		if( granger_network ){
			for( i=0; i<N; ++i) for( j=0; j<N; ++j) for( k=0; k<P; ++k ) alpha_update(alpha,i,j,k,y_aux,x_aux,rho,c,y,lambda,alpha_weights,T,N,P);
		}

		// RHO Step
		if( parcorr_network ){
			for( i=0; i<N; ++i) for( j=0; j<i; ++j ) rho_update(rho,i,j,y_aux,x_aux,alpha,c,y,lambda,rho_weights,T,N,P);
		}

		delta = 0;
		if( granger_network ) for( i=0; i<N*N*P ;    ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
		if( parcorr_network ) for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]);     }
		if( delta<toll ) break;

	}

	if( verbose ) nets_log(alpha,rho,y_aux,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,0,delta);  
}

// ALPHA update
void alpha_update(double *alpha, int i, int j, int k, double *y_aux, double *x_aux, double *rho, double *c, double **y, double lambda, double *alpha_weights, int T, int N, int P){

	int ip, jp, kp, l, t;
	double c_yx = 0;
	double c_xx = 0;
	double kappa;

	for( ip=0; ip<N; ++ip ){
		for( t=P; t<T; ++t ){
			kappa = (ip==i)? 1.0 : (-rho[ RHOIDX(ip,i) ] * sqrt(c[ip]/c[i]));
			
			y_aux[ ip*T+t ] += alpha[ ALPIDX(i,j,k,N,P) ] * kappa * y[t-k-1][j]; 
			x_aux[ ip*T+t ]  = kappa * y[t-k-1][j];

			c_yx += y_aux[ ip*T+t ] * x_aux[ ip*T+t ];
			c_xx += x_aux[ ip*T+t ] * x_aux[ ip*T+t ];
		}
	}

	alpha[ ALPIDX(i,j,k,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]);
	//Rprintf("EXACT: %d %d %d -> %f %f : beta_ls %f beta_lasso %f\n",1+i,1+j,1+k,c_yx,c_xx,c_yx/c_xx,alpha[ ALPIDX(i,j,k,N,P) ]);

	if( alpha[ ALPIDX(i,j,k,N,P) ] != 0.0 )
	{
		for( ip=0; ip<N; ++ip ){
			for( t=P; t<T; ++t ){
				y_aux[ ip*T+t ] -= alpha[ ALPIDX(i,j,k,N,P) ] * x_aux[ ip*T+t ];
			}
		}
	}
}

// ALPHA update init
void alpha_update_init(double *alpha, int i, int j, int k, double **y, double lambda, double *alpha_weights, int T, int N, int P){

	int t;
	double c_yx = 0;
	double c_xx = 0;

	for( t=P; t<T; ++t ){
		c_yx +=     y[t][i] * y[t-k-1][j];
		c_xx += y[t-k-1][j] * y[t-k-1][j];
	}

	alpha[ ALPIDX(i,j,k,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]);
}

// RHO update
void rho_update(double *rho, int i, int j, double *y_aux, double *x_aux, double *alpha, double *c, double **y, double lambda, double *rho_weights, int T, int N, int P){

	int t, ip, k;
	double c_yx = 0;
	double c_xx = 0;
	
	memset(x_aux,0,sizeof(double)*N*T);				
	for( t=P; t<T; ++t ){
		
		x_aux[ i*T+t ] = y[t][j];
		x_aux[ j*T+t ] = y[t][i];

		// resid
		for( ip=0; ip<N; ++ip ){
			for( k=0; k<P; ++k ){
				x_aux[ i*T+t ] -= alpha[ ALPIDX(j,ip,k,N,P) ] * y[t-k-1][ip];
				x_aux[ j*T+t ] -= alpha[ ALPIDX(i,ip,k,N,P) ] * y[t-k-1][ip];
			}
		}
		x_aux[ i*T+t ] = sqrt(c[j]/c[i]) * x_aux[ i*T+t ];
		x_aux[ j*T+t ] = sqrt(c[i]/c[j]) * x_aux[ j*T+t ];

		y_aux[ i*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ i*T+t ];
		y_aux[ j*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ j*T+t ];
	
		c_yx += y_aux[ i*T+t ] * x_aux[ i*T+t ] + y_aux[ j*T+t ] * x_aux[ j*T+t ];
		c_xx += x_aux[ i*T+t ] * x_aux[ i*T+t ] + x_aux[ j*T+t ] * x_aux[ j*T+t ];
	}

	// update alpha
	rho[ RHOIDX(i,j) ] = soft_thresholding(c_yx,c_xx,lambda*rho_weights[RHOIDX(i,j)] ); 

	if( rho[ RHOIDX(i,j) ] != 0.0 )
	{
		for( t=P; t<T; ++t ){          
			y_aux[ i*T+t ] -= rho[ RHOIDX(i,j) ] * x_aux[ i*T+t ];
			y_aux[ j*T+t ] -= rho[ RHOIDX(i,j) ] * x_aux[ j*T+t ];
		}
	}
}

// RHO update init
void rho_update_init(double *rho, int i, int j, double *y_aux, double *x_aux, double *alpha, double *c, double **y, double lambda, double *rho_weights, int T, int N, int P){

	int t, ip, k;
	double c_yx = 0;
	double c_xx = 0;
	
	memset(x_aux,0,sizeof(double)*N*T);				
	for( t=P; t<T; ++t ){
		//c_yx += y_aux[ i*T+t ] * x_aux[ i*T+t ] + y_aux[ j*T+t ] * x_aux[ j*T+t ];
		//c_xx += x_aux[ i*T+t ] * x_aux[ i*T+t ] + x_aux[ j*T+t ] * x_aux[ j*T+t ];
	}

	// update alpha
	//rho[ RHOIDX(i,j) ] = soft_thresholding(c_yx,c_xx,lambda*rho_weights[RHOIDX(i,j)] ); 
	//Rprintf("%d,%d > c_yx %f c_xx %f lambda %f > rho %f\n",i,j,c_yx,c_xx,lambda,rho[RHOIDX(i,j)]);
}

void nets_log(double *alpha, double *rho, double *y_aux, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int GN, int CN, int iter, double delta){

	int i;
	double rss = 0;
	double pen = 0;
	double nnz = 0;

	for( i=0; i<N*T; ++i ){
		rss += y_aux[i]*y_aux[i];
	}
      
	for( i=0; i<N*N*P; ++i ){     
        	pen += lambda * alpha_weights[i] * fabs( alpha[i] );
		nnz+= fabs(alpha[i])>0?1:0;
	}
	for( i=0; i<N*(N-1)/2; ++i ){ 
        	pen += lambda * rho_weights[i]   * fabs( rho[i] ); 
	        nnz += fabs(rho[i])>0?1:0;
	}

	if( iter > 0 ) Rprintf(" Iter: %4.4d"   , iter );  
	else Rprintf(" Converged!");  

	Rprintf(" RSS: %4.4f"    , rss/(T*N) );
	Rprintf(" Pen: %4.4f"    , pen/(T*N) );
	Rprintf(" Obj: %4.4f"    , (rss+pen)/(T*N) );
	Rprintf(" Spars: %4.4f"  , nnz / (GN*N*N*P + CN*N*(N-1)/2) );

	if( iter!= 1 ) Rprintf(" Delta: %4.4f", delta );
	Rprintf("\n");
}

void nets_sanity_check(double **y, double *alpha, double *rho, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int GN, int CN){

	int t,i,j,k;
	
	Rprintf( "NETS: T %d N %d P %d GN %d CN %d\n\n",T,N,P,GN,CN);

	// Check A
	Rprintf( "A indices\n");
	for(k=0;k<P;++k){
		for(i=0;i<N;++i){
			for(j=0;j<N;++j){
				Rprintf( "%3.3d, " , ALPIDX(i,j,k,N,P)  );
			}
			Rprintf( "\n" );	
		}
		Rprintf( "\n" );
	}

	Rprintf( "A matrix\n");
	for(k=0;k<P;++k){
		for(i=0;i<N;++i){
			for(j=0;j<N;++j){
				Rprintf( "%3.3f, " , alpha[ ALPIDX(i,j,k,N,P) ]  );
			}
			Rprintf( "\n" );	
		}
		Rprintf( "\n" );
	}
	Rprintf( "alpha vector\n");
	for( i=0;i<N*N*P; ++i) Rprintf( "%f, " , alpha[ i ]  );
	Rprintf( "\n\n" );
	Rprintf( "alpha weights\n");
	for( i=0;i<N*N*P; ++i) Rprintf( "%f, " , alpha_weights[ i ]  );
	Rprintf( "\n" );
	Rprintf( "\n" );

	// Check C
	Rprintf( "C indices\n");
	for(i=0;i<N;++i){
		for(j=0;j<N;++j){
			if( i!=j ) Rprintf( "(%d,%d) %d, " , i, j, RHOIDX(i,j)  );
			else Rprintf("(x,x) x, ");
		}
		Rprintf( "\n" );
	}
	Rprintf( "\n" );

	Rprintf( "rho vector and rho vector weights:\n" );
	for( i=0;i<N*(N-1)/2; ++i) Rprintf( "% 4.4f, " , rho[ i ]  );
	Rprintf( "\n" );
	for( i=0;i<N*(N-1)/2; ++i) Rprintf( "% 4.4f, " , rho_weights[ i ]  );
	Rprintf( "\n" );

	Rprintf( "rho matrix:\n" );
	for(i=0;i<N;++i){
		for(j=0;j<N;++j){
			if( i!=j ) Rprintf( "% 4.4f, " , rho[ RHOIDX(i,j) ]  );
			else Rprintf( "% 4.4f, ", 1.0 );
		}
		Rprintf( "\n" );	
	}
	Rprintf( "\n" );	

	Rprintf(" y matrix:\n");
	for( t=0; t<T; ++t ){
		for( i=0; i<N; ++i){
			Rprintf( "% 4.4f, " , y[t][i] );
		}
		Rprintf( "\n" );	
	}

}

/*
// ALPHA update
void alpha_update(double *alpha, int i, int j, int k, double **C_y, double *rho, double *c, double lambda, double *alpha_weights, int T, int N, int P, double **y){

	int ip, jp, kp, l; 
	double c_yx = 0;
	double c_xx = 0;
	double kappa;

	c_yx = 0.0;
	c_xx = 0.0;
	for( ip=0; ip<N; ++ip ){

		kappa = (ip==i)? 1.0 : ( -rho[ RHOIDX(i,ip) ] * sqrt(c[ip]/c[i]) );

		c_yx += kappa * C_y[ip][(1+k)*N+j];

		for( kp=0; kp<P; ++kp ){
			for( jp=0; jp<N; ++jp ){
				c_yx += -alpha[ ALPIDX(ip,jp,kp,N,P) ] * kappa * C_y[(1+kp)*N+jp][(1+k)*N+j];
				for( l=0; l<N; ++l ){
					c_yx += ( l != ip ) *  rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * alpha[ ALPIDX(l,jp,kp,N,P) ] * kappa * C_y[(1+kp)*N+jp][(1+k)*N+j];
				}

			}
		}

		for( l=0 ; l<N; ++l ){
			c_yx += ( l != ip) * -rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * kappa * C_y[l][(1+k)*N+j];
		}

		if( ip==i ) c_yx += alpha[ ALPIDX(i,j,k,N,P) ] * kappa * C_y[(1+k)*N+j][(1+k)*N+j];
		else        c_yx += -alpha[ ALPIDX(i,j,k,N,P) ] * rho[ RHOIDX(ip,i) ] * sqrt(c[ip]/c[i]) * kappa * C_y[(1+k)*N+j][(1+k)*N+j];
		c_xx += kappa * kappa * C_y[(1+k)*N+j][(1+k)*N+j];
	}

	alpha[ ALPIDX(i,j,k,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]);

	//Rprintf("QUICK: %d %d %d -> %f %f : beta_ls %f beta_lasso %f\n",1+i,1+j,1+k,c_yx,c_xx,c_yx/c_xx,soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]));
	// compute: y_aux, x_aux, c_yx, c_xx
	double *x_aux = R_alloc( N*T , double ); 
	double *y_aux = R_alloc( N*T , double ); 
	int t;
	c_yx = 0;
	c_xx = 0;
	for( ip=0; ip<N; ++ip ){
		for( t=P; t<T; ++t ){
			y_aux[ ip*T+t ] = y[t][ip];
			for( kp=0; kp<P; ++kp ){
				for( jp=0; jp<N; ++jp ){
					y_aux[ ip*T+t ] -= alpha[ ALPIDX(ip,jp,kp,N,P) ] * y[t-kp-1][jp];
					for( l=0; l<N; ++l ){
						if( l!=ip ) y_aux[ ip*T+t ] += rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * alpha[ ALPIDX(l,jp,kp,N,P) ] * y[t-kp-1][jp];
					}
				}
			}
			for( l=0 ; l<N; ++l ){
				if( l!=ip ) y_aux[ ip*T+t ] -=  rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * y[t][l];
			}

			if( ip==i ){
				y_aux[ ip*T+t ] += +alpha[ ALPIDX(i,j,k,N,P) ] * y[t-k-1][j];
				x_aux[ ip*T+t ]  = y[t-k-1][j];
			}
			else {
				y_aux[ ip*T+t ] += -alpha[ ALPIDX(i,j,k,N,P) ] * rho[ RHOIDX(ip,i) ] * sqrt(c[ip]/c[i]) * y[t-k-1][j]; 
				x_aux[ ip*T+t ]  = -rho[ RHOIDX(ip,i) ] * sqrt(c[ip]/c[i]) * y[t-k-1][j];
			}

			c_yx += y_aux[ ip*T+t ] * x_aux[ ip*T+t ];
			c_xx += x_aux[ ip*T+t ] * x_aux[ ip*T+t ];
		}
	}

	Free( x_aux );
	Free( y_aux );
	
	// update alpha
	alpha[ ALPIDX(i,j,k,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]);

	//Rprintf("EXACT: %d %d %d -> %f %f : beta_ls %f beta_lasso %f\n",1+i,1+j,1+k,c_yx,c_xx,c_yx/c_xx,soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]));
}

// RHO update
void rho_update(double *rho, int i, int j, double **C_eps, double *c, double lambda, double *rho_weights, int T, int N, int P){
  			
	int l;
	double c_yx, c_xx;

	c_yx = sqrt(c[j]/c[i]) * C_eps[i][j] + sqrt(c[i]/c[j]) * C_eps[j][i];
	c_xx = (c[j]/c[i]) * C_eps[j][j] + (c[i]/c[j]) * C_eps[i][i];

	for( l=0; l<N; ++l) {
		if( l!=i && l!=j ){
			c_yx -= rho[ RHOIDX(i,l) ] * ( sqrt(c[j]/c[i]) * sqrt(c[l]/c[i]) * C_eps[i][l] + sqrt(c[i]/c[j]) * sqrt(c[i]/c[l]) * C_eps[i][l] );
			c_yx -= rho[ RHOIDX(j,l) ] * ( sqrt(c[j]/c[i]) * sqrt(c[l]/c[j]) * C_eps[j][l] + sqrt(c[i]/c[j]) * sqrt(c[j]/c[l]) * C_eps[j][l] );
		}
	}
	
	// update rho
	rho[ RHOIDX(i,j) ] = soft_thresholding(c_yx,c_xx,lambda*rho_weights[RHOIDX(i,j)] );
}


// NETS Standard
void nets_std(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v, int *m)
{
	// variables 
	int T, N, P;
	int granger_network, parcorr_network;
	double **y, **eps;
	double *alpha_old, *rho_old;
	double delta, toll;
	double lambda;
	int iter, maxiter;
	int verbose;
	int i, j, k, kp;
	int t;
	double **C_y;
	double **C_eps;
  
	// init
	maxiter  = *m;
	toll     = 1e-4;
	verbose  = *v;
	T = *_T;
	N = *_N;
 	P = *_P;
	lambda = *_lambda;
	granger_network = *GN;
	parcorr_network = *CN;

	if( granger_network==0 ){ P = 0; };

	// allocate memory
	y   = R_alloc(T,double*);
	eps = R_alloc(T,double*);
	for( t=0; t<T; t++) { 
		y[t]   = R_alloc(N,double); 
		eps[t] = R_alloc(N,double); 
		for( i=0; i<N; i++ ) {
			y[t][i]   = _y[ T * i + t ];
			eps[t][i] = 0;
		}
	}

	if( granger_network ){
		alpha_old = R_alloc(N*N*P,double);

		C_y = R_alloc(N*(P+1),double*);
		for( k=0; k<N*(P+1); ++k){
			C_y[k] = R_alloc(N*(P+1),double);
		}

		for( k=0; k<P+1; ++k ){
			for( kp=0; kp<P+1; ++kp ){
				for( i=0; i<N; ++i ){
					for( j=0; j<N; ++j ){
						C_y[(k*N)+i][(kp*N)+j] = 0.0;
						for( t=P; t<T; ++t ) C_y[(k*N)+i][(kp*N)+j] += y[t-k][i]*y[t-kp][j];
					}
				}
			}
		}
	}
	
	if( parcorr_network ){
		rho_old   = R_alloc(N*(N-1)/2,double);
		
		C_eps = R_alloc(N,double*);
		for( i=0; i<N; ++i) C_eps[i] = R_alloc(N,double); 

		// Covariance eps
		for( i=0; i<N; ++i ){
			for( t=P; t<T; ++t ){
				eps[t][i] = y[t][i];
				for( j=0; j<N; ++j ){ for( k=0; k<P; ++k ){ eps[t][i] -= alpha[ ALPIDX(i,j,k,N,P) ] * y[t-k-1][j]; } }
			}
		}
		for( i=0; i<N; ++i ){
			for( j=0; j<=i; ++j ){
				C_eps[i][j] = 0;
				for( t=P; t<T; ++t ) C_eps[i][j] += eps[t][i]*eps[t][j];
				C_eps[j][i] = C_eps[i][j];
			}
		}
	}

	//nets_sanity_check(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,C_y,C_eps);

	// main loop
	for( iter=1; iter<=maxiter; ++iter ){

		if( verbose ) nets_log(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,iter,delta);

		// GLOBAL Update 
		if( granger_network ) memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
		if( parcorr_network ) memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

		// ALPHA Step
		if( granger_network ){
			for( i=0; i<N; ++i) for( j=0; j<N; ++j) for( k=0; k<P; ++k ) alpha_update(alpha,i,j,k,C_y,rho,c,lambda,alpha_weights,T,N,P,y);
		}
		
		if( granger_network && parcorr_network ){ 
			for( i=0; i<N; ++i ){
				for( t=P; t<T; ++t ){
					eps[t][i] = y[t][i];
					for( j=0; j<N; ++j ){ for( k=0; k<P; ++k ){ eps[t][i] -= alpha[ ALPIDX(i,j,k,N,P) ] * y[t-k-1][j]; } }
				}
			}
			for( i=0; i<N; ++i ){
				for( j=0; j<=i; ++j ){
					C_eps[i][j] = 0;
					for( t=P; t<T; ++t ) C_eps[i][j] += eps[t][i]*eps[t][j];
					C_eps[j][i] = C_eps[i][j];
				}
			}
		}

		// RHO Step
		if( parcorr_network ){
			for( i=0; i<N; ++i) for( j=0; j<i; ++j ) rho_update(rho,i,j,C_eps,c,lambda,rho_weights,T,N,P);
		}

		// Convergence Check
		delta = 0;
		if( granger_network ) for( i=0; i<N*N*P ;    ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
		if( parcorr_network ) for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]);     }
		if( delta<toll ) break;
	}
	
	if( verbose ) nets_log(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,iter,delta);
  
	// clean up	
	for(t = 0; t < T; t++){ 
		Free(y[t]); 
		Free(eps[t]);
	}
	Free(y);
	Free(eps);

	if( granger_network ){
		Free(alpha_old);

		for( k=0; k<N*(P+1); ++k ){
			Free(C_y[k]);
		}
		Free(C_y);
	}
	if( parcorr_network ) {
		Free(rho_old);

		for( i=0; i<N; ++i) { 
			Free(C_eps[i]);
		}
		Free(C_eps);
	}
}

// NETS Active Set
void nets_activeset(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v, int *m)
{
	// variables 
	int T, N, P;
	int granger_network, parcorr_network;
	double **y, **eps;
	double *alpha_old, *rho_old;
	double **alpha_set, **rho_set;
	int alpha_nact, rho_nact;
	double delta, toll;
	double lambda;
	int iter1, iter2, maxiter;
	int verbose;
	int idx, i, j, k, kp;
	int t;
	double **C_y;
	double **C_eps;
  
	// init
	maxiter  = *m;
	toll     = 1e-4;
	verbose  = *v;
	T = *_T;
	N = *_N;
 	P = *_P;
	lambda = *_lambda;
	granger_network = *GN;
	parcorr_network = *CN;

	if( granger_network==0 ){ P = 0; };

	// allocate memory
	y   = R_alloc(T,double*);
	eps = R_alloc(T,double*);
	for( t=0; t<T; t++) { 
		y[t]   = R_alloc(N,double); 
		eps[t] = R_alloc(N,double); 
		for( i=0; i<N; i++ ) {
			y[t][i]   = _y[ T * i + t ];
			eps[t][i] = 0;
		}
	}

	if( granger_network ){
		alpha_old = R_alloc(N*N*P,double);

		alpha_set = R_alloc(N*N*P,double *);
		for( idx=0; idx<N*N*P; ++idx){
			alpha_set[idx] = R_alloc(3,double); 
		}


		C_y = R_alloc(N*(P+1),double*);
		for( k=0; k<N*(P+1); ++k){
			C_y[k] = R_alloc(N*(P+1),double);
		}

		for( k=0; k<P+1; ++k ){
			for( kp=0; kp<P+1; ++kp ){
				for( i=0; i<N; ++i ){
					for( j=0; j<N; ++j ){
						C_y[(k*N)+i][(kp*N)+j] = 0.0;
						for( t=P; t<T; ++t ) C_y[(k*N)+i][(kp*N)+j] += y[t-k][i]*y[t-kp][j];
					}
				}
			}
		}
	}
	
	if( parcorr_network ){
		rho_old   = R_alloc(N*(N-1)/2,double);

		rho_set = R_alloc(N*(N-1)/2,double *);
		for( idx=0; idx<N*(N-1)/2; ++idx){
			rho_set[idx] = R_alloc(2,double); 
		}
	
		C_eps = R_alloc(N,double*);
		for( i=0; i<N; ++i) C_eps[i] = R_alloc(N,double); 

		// Covariance eps
		for( i=0; i<N; ++i ){
			for( t=P; t<T; ++t ){
				eps[t][i] = y[t][i];
				for( j=0; j<N; ++j ){ for( k=0; k<P; ++k ){ eps[t][i] -= alpha[ ALPIDX(i,j,k,N,P) ] * y[t-k-1][j]; } }
			}
		}
		for( i=0; i<N; ++i ){
			for( j=0; j<=i; ++j ){
				C_eps[i][j] = 0;
				for( t=P; t<T; ++t ) C_eps[i][j] += eps[t][i]*eps[t][j];
				C_eps[j][i] = C_eps[i][j];
			}
		}
	}

	// check! 
	//nets_sanity_check(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,C_y,C_eps);

	// active set
	for( iter1=1; iter1<=maxiter; ++iter1 ){
			
		if( verbose ) nets_log(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,iter1,delta);

		// find active set
		if( granger_network ){
			alpha_nact = 0;
			for( i=0; i<N; ++i){ 
				for( j=0; j<N; ++j){ 
					for( k=0; k<P; ++k ){
						if( alpha[ ALPIDX(i,j,k,N,P) ]!=0 ){
							alpha_set[ alpha_nact ][0] = i;
							alpha_set[ alpha_nact ][1] = j;
							alpha_set[ alpha_nact ][2] = k;
							++alpha_nact;
						}
					}
				}
			}
		}
		if( parcorr_network ){
			rho_nact = 0;
			for( i=0; i<N; ++i){ 
				for( j=0; j<i; ++j){ 
					if( rho[ RHOIDX(i,j) ]!=0 ){
						rho_set[ rho_nact ][0] = i;
						rho_set[ rho_nact ][1] = j;
						++rho_nact;
					}
				}
			}
		}

		// ACTIVE Update 
		for( iter2=1; iter2<=maxiter; ++iter2 ){


			if( granger_network ) memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
			if( parcorr_network ) memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

			// ALPHA Step
			if( granger_network ){
				for( idx=0; idx<alpha_nact; ++idx){
					 i = alpha_set[idx][0];
					 j = alpha_set[idx][1];
					 k = alpha_set[idx][2];
					 alpha_update(alpha,i,j,k,C_y,rho,c,lambda,alpha_weights,T,N,P,y);
				}
			}
			
			if( granger_network && parcorr_network ){ 
				for( i=0; i<N; ++i ){
					for( t=P; t<T; ++t ){
						eps[t][i] = y[t][i];
						for( j=0; j<N; ++j ){ for( k=0; k<P; ++k ){ eps[t][i] -= alpha[ ALPIDX(i,j,k,N,P) ] * y[t-k-1][j]; } }
					}
				}
				for( i=0; i<N; ++i ){
					for( j=0; j<=i; ++j ){
						C_eps[i][j] = 0;
						for( t=P; t<T; ++t ) C_eps[i][j] += eps[t][i]*eps[t][j];
						C_eps[j][i] = C_eps[i][j];
					}
				}
			}

			// RHO Step
			if( parcorr_network ){
				for( idx=0; idx<rho_nact; ++idx){
					i = rho_set[idx][0];
					j = rho_set[idx][1];
					rho_update(rho,i,j,C_eps,c,lambda,rho_weights,T,N,P);
				}
			}

			// Convergence Check
			delta = 0;
			if( granger_network ) for( i=0; i<N*N*P ;    ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
			if( parcorr_network ) for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]);     }
			if( delta<toll ) break;
			if( iter1==1 ) break;
		}

		// GLOBAL UPDATE
		if( granger_network ) memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
		if( parcorr_network ) memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

		// ALPHA Step
		if( granger_network ){
			for( i=0; i<N; ++i) for( j=0; j<N; ++j) for( k=0; k<P; ++k ) alpha_update(alpha,i,j,k,C_y,rho,c,lambda,alpha_weights,T,N,P,y);
		}
			
		if( granger_network && parcorr_network ){ 
			for( i=0; i<N; ++i ){
				for( t=P; t<T; ++t ){
					eps[t][i] = y[t][i];
					for( j=0; j<N; ++j ){ for( k=0; k<P; ++k ){ eps[t][i] -= alpha[ ALPIDX(i,j,k,N,P) ] * y[t-k-1][j]; } }
				}
			}
			for( i=0; i<N; ++i ){
				for( j=0; j<=i; ++j ){
					C_eps[i][j] = 0;
					for( t=P; t<T; ++t ) C_eps[i][j] += eps[t][i]*eps[t][j];
					C_eps[j][i] = C_eps[i][j];
				}
			}
		}

		// RHO Step
		if( parcorr_network ){
			for( i=0; i<N; ++i) for( j=0; j<i; ++j ) rho_update(rho,i,j,C_eps,c,lambda,rho_weights,T,N,P);
		}

		delta = 0;
		if( granger_network ) for( i=0; i<N*N*P ;    ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
		if( parcorr_network ) for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]);     }
		if( delta<toll ) break;

	}

	if( verbose ) nets_log(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,iter1,delta);
  
	// clean up	
	for(t = 0; t < T; t++){ 
		Free(y[t]); 
		Free(eps[t]);
	}
	Free(y);
	Free(eps);

	if( granger_network ){
		Free(alpha_old);

		for( idx=0;idx<N*N*P; ++idx){
			Free(alpha_set[idx]);
		}
		Free(alpha_set);

		for( k=0; k<N*(P+1); ++k ){
			Free(C_y[k]);
		}
		Free(C_y);
	}
	if( parcorr_network ) {
		Free(rho_old);

		for( idx=0;idx<N*(N-1)/2; ++idx){
			Free(rho_set[idx]);
		}
		Free(rho_set);

		for( i=0; i<N; ++i) { 
			Free(C_eps[i]);
		}
		Free(C_eps);
	}
}

void nets_sanity_check(double **y, double *alpha, double *rho, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int GN, int CN, double ***C_y, double C_eps){

	int i,j,k;
		
	Rprintf( "NETS: T %d N %d P %d GN %d CN %d\n\n",T,N,P,GN,CN);

	// Check A
	if( GN ){
		Rprintf( "A indices\n");
		for(k=0;k<P;++k){
			for(i=0;i<N;++i){
				for(j=0;j<N;++j){
					Rprintf( "%d, " , ALPIDX(i,j,k,N,P)  );
				}
				Rprintf( "\n" );	
			}
			Rprintf( "\n" );
		}
		Rprintf( "\n" );

		Rprintf( "A\n");
		for(k=0;k<P;++k){
			for(i=0;i<N;++i){
				for(j=0;j<N;++j){
					Rprintf( "%f, " , alpha[ ALPIDX(i,j,k,N,P) ]  );
				}
				Rprintf( "\n" );	
			}
			Rprintf( "\n" );
		}
		for( i=0;i<N*N*P; ++i) Rprintf( "%f, " , alpha[ i ]  );
		Rprintf( "\n" );
	}
	Rprintf( "\n" );

	// Check R
	if( CN ){
		Rprintf( "C indices\n");
		for(i=0;i<N;++i){
			for(j=0;j<N;++j){
				if( i!=j ) Rprintf( "%d, " , RHOIDX(i,j)  );
				else Rprintf("x, ");
			}
			Rprintf( "\n" );
		}
		Rprintf( "\n" );

		Rprintf( "rho vector and rho vector weights:\n" );
		for( i=0;i<N*(N-1)/2; ++i) Rprintf( "% 4.4f, " , rho[ i ]  );
		Rprintf( "\n" );
		for( i=0;i<N*(N-1)/2; ++i) Rprintf( "% 4.4f, " , rho_weights[ i ]  );
		Rprintf( "\n" );

		Rprintf( "rho matrix:\n" );
		for(i=0;i<N;++i){
			for(j=0;j<N;++j){
				if( i!=j ) Rprintf( "% 4.4f, " , rho[ RHOIDX(i,j) ]  );
				else Rprintf( "% 4.4f, ", 1.0 );
			}
			Rprintf( "\n" );	
		}
		Rprintf( "\n" );	
	}

	// check C_y
	if( GN ){
		for(k=0; k<P+1; ++k){
			Rprintf( "C(%d)\n" , k );
			for(i=0;i<N;++i){
				for(j=0;j<N;++j){
					Rprintf( "% 4.4f, " , C_y[k][i][j]  );
				}
				Rprintf( "\n" );	
			}	
		}
		Rprintf( "\n" );
	}
}

void nets_log(double **y, double *alpha, double *rho, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int GN, int CN, int iter, double delta){

	int i;
	double pen = 0;
	double nnz = 0;
      
	for( i=0; i<N*N*P; ++i ){     
        	pen += lambda * alpha_weights[i] * fabs( alpha[i] );
		nnz+= fabs(alpha[i])>0?1:0;
	}
	for( i=0; i<N*(N-1)/2; ++i ){ 
        	pen += lambda * rho_weights[i]   * fabs( rho[i] ); 
	        nnz += fabs(rho[i])>0?1:0;
	}

	Rprintf(" Iter: %4.4d"   , iter );  
	Rprintf(" RSS: %4.4f"    , 0.0/(T*N) );
	Rprintf(" Pen: %4.4f"    , pen/(T*N) );
	Rprintf(" Obj: %4.4f"    , 0.0/(T*N) );
	Rprintf(" Spars: %4.4f"  , nnz / (GN*N*N*P + CN*N*(N-1)/2) );
	Rprintf(" Delta: %4.4f\n", delta );
}

*/
