
#include <R.h>
#include <math.h>

#define MAX(a,b)          (a>b?a:b)
#define MIN(a,b)          (a<b?a:b)
#define RHOIDX(i,j)       ( MAX(i,j)*(MAX(i,j)-1)/2 + MIN(i,j) )
#define ALPIDX(i,j,p,N,P) ( p*N*N + i*N + j )

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

// ALPHA update
void alpha_update(double *alpha, int i, int j, int k, double ***C_y, double *rho, double *c, double lambda, double *alpha_weights, int T, int N, int P, double **y){

	int ip, jp, kp, l; 
	double c_yx = 0;
	double c_xx = 0;
	double c_tmp;

	c_yx = 0.0;
	c_xx = 0.0;
	for( ip=0; ip<N; ++ip ){

		if( ip == i ) c_yx += C_y[1+k][i][j];
		else          c_yx += -rho[ RHOIDX(i,ip) ] * sqrt(c[ip]/c[i]) * C_y[1+k][ip][j];

		for( kp=0; kp<P; ++kp ){
			for( jp=0; jp<N; ++jp ){

				c_tmp = ( (k-kp)>=0 ) ? ( C_y[k-kp][jp][j] ):( C_y[kp-k][j][jp] );

				if( ip == i ){
					c_yx += -alpha[ ALPIDX(ip,jp,kp,N,P) ] * c_tmp * ( (double) ( ip!=i || jp!=j || kp!=k ) );
				}
				else {
					c_yx += alpha[ ALPIDX(ip,jp,kp,N,P) ] * rho[ RHOIDX(i,ip) ] * sqrt(c[ip]/c[i])  * c_tmp;
				}

				for( l=0; l<N; ++l ){
					if( l!=ip ){
						if( ip == i ) c_yx +=  
						rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * alpha[ ALPIDX(l,jp,kp,N,P) ] * c_tmp * ( (double) ( l!=i || jp!=j || kp!=k ) );
						else          c_yx += 
						-rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * alpha[ ALPIDX(l,jp,kp,N,P) ] * rho[ RHOIDX(i,ip) ] * sqrt(c[ip]/c[i]) * c_tmp * ( (double) ( l!=i || jp!=j || kp!=k ) );
					}
				}

			}
		}
		for( l=0 ; l<N; ++l ){
			if( l!=ip ){
				if( ip==i ) c_yx += -rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * C_y[1+k][l][j];
				else        c_yx +=  rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * rho[ RHOIDX(i,ip) ] * sqrt(c[ip]/c[i]) * C_y[1+k][l][j];
			}
		}

		if( ip==i ){
			//c_yx += +alpha[ ALPIDX(i,j,k,N,P) ] * C_y[0][j][j];
			c_xx += C_y[0][j][j];
		}
		else {
			c_xx += rho[ RHOIDX(ip,i) ] * rho[ RHOIDX(ip,i) ] * (c[ip]/c[i]) * C_y[0][j][j];
		}
	}
	
	// update alpha
	//Rprintf("%d %d %d -> %f %f : beta_ls %f beta_lasso %f\n",1+i,1+j,1+k,c_yx,c_xx,c_yx/c_xx,soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]));
	//alpha[ ALPIDX(i,j,k,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]);

	// compute: y_aux, x_aux, c_yx, c_xx
	/*
	double *x_aux = Calloc( N*T , double ); 
	double *y_aux = Calloc( N*T , double ); 
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
						if( l != ip ) y_aux[ ip*T+t ] += rho[ RHOIDX(ip,l) ] * sqrt(c[l]/c[ip]) * alpha[ ALPIDX(l,jp,kp,N,P) ] * y[t-kp-1][jp];
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
	*/

	// update alpha
	alpha[ ALPIDX(i,j,k,N,P) ] = soft_thresholding(c_yx,c_xx,lambda*alpha_weights[ALPIDX(i,j,k,N,P)]);

	//Rprintf("%d %d %d -> %f %f : beta_ls %f beta_lasso %f\n",1+i,1+j,1+k,c_yx,c_xx,c_yx/c_xx,alpha[ ALPIDX(i,j,k,N,P) ]);

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


//
void nets_std(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v)
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
	double ***C_y;
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

		C_y = Calloc(N,double**);
		for( k=0; k<P+1; ++k){
			C_y[k] = Calloc(N,double*);
			for( i=0; i<N; ++i) C_y[k][i] = Calloc(N,double); 
		}

		// Covariance y
		for( k=0; k<P+1; ++k ){
			for( i=0; i<N; ++i ){
				for( j=0; j<N; ++j ){
					C_y[k][i][j] = 0;
					for( t=P; t<T; ++t ) C_y[k][i][j] += y[t][i]*y[t-k][j];
				}
			}
		}
	}
	
	if( parcorr_network ){
		rho_old   = Calloc(N*(N-1)/2,double);
		
		C_eps = Calloc(N,double*);
		for( i=0; i<N; ++i) C_eps[i] = Calloc(N,double); 

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
	nets_sanity_check(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,granger_network,parcorr_network,C_y,C_eps);

	// main loop
	delta = 0.0;
	for( iter=1; iter<=maxiter; ++iter ){

		if( verbose ) nets_log(y,alpha,rho,lambda,alpha_weights,rho_weights,T,N,P,iter,delta);

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
  
	// clean up	
	for(t = 0; t < T; t++){ 
		Free(y[t]); 
		Free(eps[t]);
	}
	Free(y);
	Free(eps);

	if( granger_network ){
		Free(alpha_old);

		for( k=0; k<P+1; ++k) {
			for( i=0; i<N; ++i) { 
				Free(C_y[k][i]);
			}
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

// NETS - LONG FORM
void nets_long_log(double *alpha, double *rho, double *y_aux, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int iter, double delta){

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

void nets_long(double *alpha, double *rho, double *alpha_weights, double *rho_weights, double *_lambda, double *_y, int *_T, int *_N, int *_P, double *c, int *GN, int *CN, int *v) {

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
	toll     = 1e-6;
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
		      nets_long_log(alpha,rho,y_aux,lambda,alpha_weights,rho_weights,T,N,P,iter,delta);
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

void nets_sanity_check(double **y, double *alpha, double *rho, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int GN, int CN, double ***C_y, double C_eps){

	int i,j,k;
		
	Rprintf( "NETS: T %d N %d P %d GN %d CN %d\n\n",T,N,P,GN,CN);

	// Check A
	if( GN ){
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

void nets_log(double **y, double *alpha, double *rho, double lambda, double *alpha_weights, double *rho_weights, int T, int N, int P, int iter, double delta){

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
	Rprintf(" Spars: %4.4f"  , nnz/(N*N*P+N*(N-1)/2) );
	Rprintf(" Delta: %4.4f\n", delta );
}

