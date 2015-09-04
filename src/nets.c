
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

void alpha_update(double *alpha, int i, int j, int p, double *y_aux, double *x_aux, double **y, double *rho, double *c, double lambda, double *alpha_weights, int T, int N, int P){

	int k, t;  				
	double c_yx = 0;
	double c_xx = 0;

	for( k=0; k<N; ++k ){
		for( t=P; t<T; ++t ){
			if( k==i ){
				x_aux[ k*T+t ] = y[t-p-1][j];
			}
			else{
				x_aux[ k*T+t ] = -rho[ RHOIDX(i,k) ] * sqrt(c[k]/c[i]) * y[t-p-1][j];
			}

			y_aux[ k*T+t ] += alpha[ ALPIDX(i,j,p,N,P) ]*x_aux[ k*T+t ];

			c_yx -= y_aux[ k*T+t ] * x_aux[ k*T+t ];
			c_xx += x_aux[ k*T+t ] * x_aux[ k*T+t ];
		}
	}

	// update alpha
	alpha[ ALPIDX(i,j,p,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,p,N,P)]);

	// update y_aux (if new coeff is != 0) 
	if( alpha[ ALPIDX(i,j,p,N,P) ] != 0.0 ){
		for( k=0; k<N; ++k ){
			for( t=P; t<T; ++t ){          
				y_aux[ k*T+t ] -= alpha[ ALPIDX(i,j,p,N,P) ]*x_aux[ k*T+t ];
			}
		}
	}
}

void rho_update(double *rho, int i, int j, double *y_aux, double *x_aux, double **eps, double *c, double lambda, double *rho_weights, int T, int N, int P){
  			
	int t;
	double c_yx, c_xx;

	// compute: y_aux, x_aux, c_yx, c_xx 
	c_yx = 0;
	c_xx = 0;
	
	memset(x_aux,0,sizeof(double)*N*T);				
	for( t=P; t<T; ++t ){
		
		x_aux[ i*T+t ] = sqrt(c[j]/c[i]) * eps[t][j];
		x_aux[ j*T+t ] = sqrt(c[i]/c[j]) * eps[t][i];

		y_aux[ i*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ i*T+t ];
		y_aux[ j*T+t ] += rho[ RHOIDX(i,j) ] * x_aux[ j*T+t ];
	
		c_yx -= y_aux[ i*T+t ] * x_aux[ i*T+t ] + y_aux[ j*T+t ] * x_aux[ j*T+t ];
		c_xx += x_aux[ i*T+t ] * x_aux[ i*T+t ] + x_aux[ j*T+t ] * x_aux[ j*T+t ];
	}

	Rprintf("%f %f\n",c_yx,c_xx);

	// update rho
	rho[ RHOIDX(i,j) ] = soft_thresholding(c_yx,c_xx,lambda*rho_weights[RHOIDX(i,j)] ); 

	// update y_aux (if new coeff is != 0) 
	if( rho[ RHOIDX(i,j) ] != 0.0 ){
		for( t=P; t<T; ++t ){          
			y_aux[ i*T+t ] -= rho[ RHOIDX(i,j) ] * x_aux[ i*T+t ];
			y_aux[ j*T+t ] -= rho[ RHOIDX(i,j) ] * x_aux[ j*T+t ];
		}
	}
}

void nets_log(double *alpha, double *rho, double *y_aux, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int iter, double delta){

	int i;
	double rss = 0;
	double pen = 0;
	int nnz= 0;
      
	for( i=0; i<N*T; ++i ) rss += (y_aux[i])*(y_aux[i]);
	for( i=0; i<N*N*P; ++i ){     
        	pen += lambda * alpha_weights[i] * fabs( alpha[i] );
		nnz+= fabs(alpha[i])>0?1:0;
	}
	for( i=0; i<N*(N-1)/2; ++i ){ 
        	pen += lambda * rho_weights[i]   * fabs( rho[i] ); 
	        nnz+= fabs(rho[i])>0?1:0;
	}

	Rprintf("Iter: %3.3d Rss %3.3f Pen %3.3f Obj %3.3f Spars %d/%d  Delta: %f\n",iter,rss/(T*N),pen/(T*N),(rss+pen)/(T*N),nnz,N*N*P+N*(N-1)/2,delta);  
}


//
void nets2(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v)
{
	// variables 
	int T, N, P;
	double **y, **eps;
	double *y_aux, *x_aux;
	double *alpha_old, *rho_old;
	double **alpha_activeset, **rho_activeset;
	int naa, nra;
	double delta, toll;
	double lambda;
	int iter1, iter2, maxiter;
	int verbose;
	int t, i, j, k, p, idx;
	int granger_network;
	int parcorr_network;
  
	// init
	maxiter  = 50;
	toll     = 1e-4;
	verbose  = *v;
	T = *_T;
	N = *_N;
 	P = *_P;
	lambda = *_lambda;
	granger_network = *GN;
	parcorr_network = *CN;

	Rprintf("nets v2\n");  

	// allocate memory
	y   = Calloc(T,double*);
	eps = Calloc(T,double*);
	for( t=0; t<T; t++) { 
		y[t]   = Calloc(N,double); 
		eps[t] = Calloc(N,double); 
		for( i=0; i<N; i++ ) {
			y[t][i]   = _y[ T * i + t ];
			eps[t][i] = 0;
		}
	}

	if( granger_network ){
		alpha_old = Calloc(N*N*P,double);

		naa = 0;
		alpha_activeset = Calloc(N*N*P,double*);
		for( i=0; i<N*N*P; i++) { 
			alpha_activeset[i] = Calloc(3,double); 
			alpha_activeset[i][0] = 0;
			alpha_activeset[i][1] = 0;
			alpha_activeset[i][2] = 0;
		}
	}
	
	if( parcorr_network ){
		rho_old   = Calloc(N*(N-1)/2,double);

		nra = 0;
		rho_activeset = Calloc(N*(N-1)/2,double*);
		for( i=0; i<N*(N-1)/2; i++) { 
			rho_activeset[i] = Calloc(2,double); 
			rho_activeset[i][0] = 0;
			rho_activeset[i][1] = 0;
		}
	}

	y_aux     = Calloc(T*N,double);
	x_aux     = Calloc(T*N,double);

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
					y_aux[i*T+t] -= rho[ RHOIDX(i,j) ] * sqrt(c[i]/c[j]) * y[t][j];
					for( p=0; p<P; ++p ){
						for( k=0; k<N; ++k ){
							y_aux[i*T+t] += rho[ RHOIDX(i,j) ] * sqrt(c[i]/c[j]) * alpha[ ALPIDX(j,k,p,N,P) ] * y[t-p-1][k];
						}
					}
				}
			}
		}
	}

	// main loop
	delta = 0.0;
	for( iter1=1; iter1<=maxiter; ++iter1 ){

		if( verbose ) nets_log(alpha,rho,y_aux,lambda,alpha_weights,rho_weights,T,N,P,iter1,delta);

		// ACTIVE Set Update
		naa = 0;	
		for( i=0; i<N; ++i) for( j=0; j<N; ++j) for( p=0; p<P; ++p ){
			if( fabs(alpha[ALPIDX(i,j,p,N,P)]) > 0.0 ){
				alpha_activeset[naa][0] = i;
				alpha_activeset[naa][1] = j;
				alpha_activeset[naa][2] = p;
				naa++;
			}
		}

		nra = 0;	
		for( i=0; i<N; ++i) for( j=0; j<i; ++j) {
			if( fabs(rho[RHOIDX(i,j)]) > 0.0 ){
				rho_activeset[nra][0] = i;
				rho_activeset[nra][1] = j;
				nra++;
			}
		}

		for( iter2=1; iter2<=10; ++iter2 ){
			memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
			memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

			for( idx=0; idx<naa; ++idx){ 
				i = alpha_activeset[idx][0];
				j = alpha_activeset[idx][1];
				p = alpha_activeset[idx][2];
				alpha_update(alpha,i,j,p,y_aux,x_aux,y,rho,c,lambda,alpha_weights,T,N,P);
			}
			for( i=0; i<N; ++i ){
				for( t=P; t<T; ++t ){
					eps[t][i] = y[t][i];
					for( j=0; j<N; ++j ){ for( p=0; p<P; ++p ){ eps[t][i] -= alpha[ ALPIDX(i,j,p,N,P) ] * y[t-p-1][j]; } }
				}
			}

			for( idx=0; idx<nra; ++idx){ 
				i = rho_activeset[idx][0];
				j = rho_activeset[idx][1];
				rho_update(rho,i,j,y_aux,x_aux,eps,c,lambda,rho_weights,T,N,P);
			}

			delta = 0;
			for( i=0; i<N*N*P ; ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
			for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]); }
			if( delta<toll ) break;
		}

		// GLOBAL Update 
		memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
		memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

		// ALPHA Step
		for( i=0; i<N; ++i) for( j=0; j<N; ++j) for( p=0; p<P; ++p ) alpha_update(alpha,i,j,p,y_aux,x_aux,y,rho,c,lambda,alpha_weights,T,N,P);
		
		// RHO Step
		for( i=0; i<N; ++i ){
			for( t=P; t<T; ++t ){
				eps[t][i] = y[t][i];
				for( j=0; j<N; ++j ){ for( p=0; p<P; ++p ){ eps[t][i] -= alpha[ ALPIDX(i,j,p,N,P) ] * y[t-p-1][j]; } }
			}
		}
    
		for( i=0; i<N; ++i) for( j=0; j<i; ++j ) rho_update(rho,i,j,y_aux,x_aux,eps,c,lambda,rho_weights,T,N,P);

		// Convergence Check
		delta = 0;
		for( i=0; i<N*N*P ; ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
		for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]); }
		if( delta<toll ) break;
	}
  
	// clean up	
	for(t = 0; t < T; t++){ 
		Free(y[t]); 
		Free(eps[t]);
	}
	Free(y);
	Free(eps);

	if( granger_network ){
		Free(alpha_old);
		for(i = 0; i < N*N*P; i++) Free(alpha_activeset[i]); 
		Free(alpha_activeset);
	}
	if( parcorr_network ){
		Free(rho_old);
		for(i = 0; i < N*(N-1)/2; i++) Free(rho_activeset[i]); 
		Free(rho_activeset);
	}

	Free(y_aux);
	Free(x_aux);
}


//
void nets(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v) {

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
	double c_yx, c_xx;
  
	// init
	maxiter  = 100;
	toll     = 1e-4;
	verbose  = *v;
	T = *_T;
	N = *_N;
 	P = *_P;
	lambda = *_lambda;

	Rprintf("nets v1\n");  

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
					y_aux[i*T+t] -= rho[ RHOIDX(i,j) ] * sqrt(c[i]/c[j]) * y[t][j];
					for( p=0; p<P; ++p ){
						for( k=0; k<N; ++k ){
							y_aux[i*T+t] += rho[ RHOIDX(i,j) ] * sqrt(c[i]/c[j]) * alpha[ ALPIDX(j,k,p,N,P) ] * y[t-p-1][k];
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
      nets_log(alpha,rho,y_aux,lambda,alpha_weights,rho_weights,T,N,P,iter,delta);
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
								x_aux[ k*T+t ] = -rho[ RHOIDX(i,k) ] * sqrt(c[k]/c[i]) * y[t-p-1][j];
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
					x_aux[ i*T+t ] = sqrt(c[j]/c[i]) * x_aux[ i*T+t ];
					x_aux[ j*T+t ] = sqrt(c[i]/c[j]) * x_aux[ j*T+t ];

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



/////////////////////
void alpha_update2(double *alpha, int i, int j, int k, double **C_y, double *rho, double *c, double lambda, double *alpha_weights, int T, int N, int P){

	int l, h, m; 
	double beta;
	double c_yx = 0;
	double c_xx = 0;

	c_yx = C_y[i][(k+1)*N+j];
	c_xx = 0;
	for( l=0; l<N; ++l ){

		// cov yx
		for(h=0; h<P; ++h){
			beta = ( l!=j && h!=k )?ALPIDX(i,l,h,N,P):0.0;
			for(m=0; m<N; ++m){
				beta -= rho[ RHOIDX(i,m) ] * sqrt(c[m]/c[i]) * ( m!=i && l!=j )?ALPIDX(m,l,k,N,P):0.0;
			}
			if( l==i ){
				c_yx -= beta * C_y[i][(k+1)*P+j];
			}
			else{
				c_yx -= beta * (-rho[ RHOIDX(i,l) ] * sqrt(c[l]/c[i]) ) * C_y[i][(k+1)*P+j];
			}
		}	

		if( l != i ) c_yx -= rho[ RHOIDX(i,l) ]*rho[ RHOIDX(i,l) ] * c[l]/c[i] * C_y[i][l];

		// cov xx
		if( l==i ){
			c_xx += C_y[j][j];
		}
		else{
			c_xx += rho[ RHOIDX(i,l) ]*rho[ RHOIDX(i,l) ] * c[l]/c[i] * C_y[j][j];
		}
	}

	c_yx = -c_yx;
	//Rprintf("%f %f\n",c_yx,c_xx);

	// update alpha
	alpha[ ALPIDX(i,j,k,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]);
}



void rho_update2(double *rho, int i, int j, double **C_eps, double *c, double lambda, double *rho_weights, int T, int N, int P){
  			
	int k, t;
	double c_yx, c_xx;

	c_yx = sqrt(c[j]/c[i]) * C_eps[i][j] + sqrt(c[i]/c[j]) * C_eps[j][i];
	c_xx = (c[j]/c[i]) * C_eps[j][j] + (c[i]/c[j]) * C_eps[i][i];

	for( k=0; k<N; ++k) {
		if( k!=i && k!=j ){
			c_yx -= rho[ RHOIDX(i,k) ] * ( sqrt(c[j]/c[i]) * sqrt(c[k]/c[i]) * C_eps[i][k] + sqrt(c[i]/c[j]) * sqrt(c[i]/c[k]) * C_eps[i][k] );
			c_yx -= rho[ RHOIDX(j,k) ] * ( sqrt(c[j]/c[i]) * sqrt(c[k]/c[j]) * C_eps[j][k] + sqrt(c[i]/c[j]) * sqrt(c[j]/c[k]) * C_eps[j][k] );
		}
	}

	c_yx = -c_yx;
	
	// update rho
	rho[ RHOIDX(i,j) ] = soft_thresholding(c_yx,c_xx,lambda*rho_weights[RHOIDX(i,j)] );
}

void nets_log2(double *alpha, double *rho, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int iter, double delta){

	int i;
	double rss = 0;
	double pen = 0;
	int nnz= 0;
      
	for( i=0; i<N*N*P; ++i ){     
        	pen += lambda * alpha_weights[i] * fabs( alpha[i] );
		nnz+= fabs(alpha[i])>0?1:0;
	}
	for( i=0; i<N*(N-1)/2; ++i ){ 
        	pen += lambda * rho_weights[i]   * fabs( rho[i] ); 
	        nnz+= fabs(rho[i])>0?1:0;
	}

	Rprintf("Iter: %3.3d Rss %3.3f Pen %3.3f Obj %3.3f Spars %d/%d  Delta: %f\n",iter,1,pen/(T*N),(1+pen)/(T*N),nnz,N*N*P+N*(N-1)/2,delta);  
}

//
void nets3(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v)
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
	int i, j, k;
	int t;
	double **C_y;
	double **C_eps;
  
	// init
	maxiter  = 100;
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
	y   = Calloc(T,double*);
	eps = Calloc(T,double*);
	for( t=0; t<T; t++) { 
		y[t]   = Calloc(N,double); 
		eps[t] = Calloc(N,double); 
		for( i=0; i<N; i++ ) {
			y[t][i]   = _y[ T * i + t ];
			eps[t][i] = 0;
		}
	}

	if( granger_network ){
		alpha_old = Calloc(N*N*P,double);

		C_y = Calloc(N,double*);
		for( i=0; i<N; ++i) C_y[i] = Calloc(N*(P+1),double); 
	}
	
	if( parcorr_network ){
		rho_old   = Calloc(N*(N-1)/2,double);
		
		C_eps = Calloc(N,double*);
		for( i=0; i<N; ++i) C_eps[i] = Calloc(N,double); 
	}

	// Covariance y 
	for( i=0; i<N; ++i ){
		for( j=0; j<N; ++j ){
			for( k=0; k<P+1; ++k ){
				C_y[i][k*N+j] = 0;
				for( t=P; t<T; ++t ) C_y[i][k*N+j] += y[t][i]*y[t-k][j];
			}
		}
	}

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

	// main loop
	delta = 0.0;
	for( iter=1; iter<=maxiter; ++iter ){

		if( verbose ) nets_log2(alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,iter,delta);

		// GLOBAL Update 
		if( granger_network ) memcpy(alpha_old,alpha,sizeof(double)*(N*N*P));
		if( parcorr_network ) memcpy(rho_old  ,rho  ,sizeof(double)*(N*(N-1)/2));

		// ALPHA Step
		if( granger_network ){

			for( i=0; i<N; ++i) for( j=0; j<N; ++j) for( k=0; k<P; ++k ) 
				alpha_update2(alpha,i,j,k,C_y,rho,c,lambda,alpha_weights,T,N,P);

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

			for( i=0; i<N; ++i) for( j=0; j<i; ++j ) 
				rho_update2(rho,i,j,C_eps,c,lambda,rho_weights,T,N,P);
		}

		// Convergence Check
		delta = 0;
		if( granger_network ) for( i=0; i<N*N*P ; ++i ){ delta += fabs(alpha[i]-alpha_old[i]); }
		if( parcorr_network ) for( i=0; i<N*(N-1)/2; ++i ){ delta += fabs(rho[i]-rho_old[i]); }
		if( delta<toll ) break;
	}
  
	// clean up	
	for(t = 0; t < T; t++){ 
		Free(y[t]); 
		Free(eps[t]);
	}
	Free(y);
	Free(eps);

	if( granger_network ){
		Free(alpha_old);

		for( i=0; i<N; ++i) { 
			Free(C_y[i]);
		}
		Free(C_y);
	}
	if( parcorr_network ){
		Free(rho_old);

		for( i=0; i<N; ++i) { 
			Free(C_eps[i]);
		}
		Free(C_eps);
	}
}


