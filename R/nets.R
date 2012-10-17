
.packageName <- "nets"

.First.lib <- function(lib, pkg){ library.dynam("nets", pkg, lib) }

.onLoad.lib <- function(lib, pkg){ library.dynam("nets", pkg, lib) }

nets <- function(){ }

plot.nets <- function(network){ }

nets.alasso <- function(y,X,lambda,w='adaptive'){
	
	M <- nrow(y)
	N <- ncol(X)

	toll <- 1e-6
	maxiter <- 20

	# check inputs
	if( any( !is.finite(y) ) ){ stop('The response vector contains non finite values.') }
	if( any( !is.finite(X) ) ){ stop('The data matrix contains non finite values.') }
	if( lambda < 0 ){ stop('The ALASSO penalty is negative') }

	# adaptive lasso weights
	if( w=='adaptive' ) {
		if( ncol(X) < nrow(y) ){
			beta.pre <- coef( lm( y ~ 0+X ) )
			w <- 1/abs(beta.pre)
		}
		else {
			# TODO: ridge
			w <- rep(1,N)
		}
	}
	else {
		w <- rep(1,N);
	}

	# call shooting algorithm
        results <- .C("shooting", theta=as.double(rep(0,N)) , as.double(y) , as.double(X) , as.double(lambda) , as.double(w) , as.integer(M) , as.integer(N) )

	# packaging results
	results <- list( theta=results$theta , eps=(y-X%*%results$theta) )
}

nets.sparse.varfit <- function(y,p=1,lambda=0,lambda.range=NULL,v=NULL,w='adaptive',verbose=FALSE){

	T <- nrow(y)
	N <- ncol(y)

	# lambda range
	if( is.null(lambda.range) )
	{
		lambda.range <- lambda
		bic <- rep(0,1)
	}
	else
	{
		lambda.range <- rev( lambda.range )
		bic <- rep( 0 , length(lambda.range) )
	}

	# THIS IS A BIT DUMB!
	cat('Sparse VAR estimation: ')
	stuff = list()
	trace = matrix(0,length(lambda.range),N*N*p)
	for( li in 1:length(lambda.range) )
	{
		cat('.')

		bic[li] <- 0
		A   <- array( 0 , dim=c(p,N,N) )
		eps <- matrix( 0 , T , N )
		for( i in 1:N )
		{
			Y <- matrix( 0 , T-p , 1 ) 
			X <- matrix( 0 , T-p , p*N )
			
			Y[] <- y[ (p+1):T , i]
			for( l in 1:p )
			{
				X[, ( (p-1)*N + p ):( p*N ) ] <- y[ (p+1-l):(T-l) , ]
			}	
		
			results <- nets.alasso( Y , X , v=v, w=w, lambda=lambda.range[li] )

			for( l in 1:p )
			{
				A[l,i,] <- results$theta[ ( (p-1)*N + p ):( p*N ) ]
			}

			eps[ (p+1):T ,i] <- results$eps

			bic[li] <- bic[li] + T * log( sum(eps^2) ) + log(T) * sum( results$theta != 0 ); 
		}

		trace[li,] <- A[1:(N*N*p)]

		stuff[[li]] <- list()
		stuff[[li]]$A <- A
		stuff[[li]]$eps <- eps
	}
	cat('\n')

	opt <- (1:length(lambda.range))[ rank(bic,ties.method="first")==1 ]
	lambda.hat <- lambda.range[ opt ]
	A <- stuff[[ opt ]]$A
	eps <- stuff[[ opt ]]$eps

	for( li in 1:length(lambda.range) )
	{
		trace[li,] <- trace[li,] / (trace[length(lambda.range),]+0.01)
	}

	list( A=A , eps=eps , lambda=lambda.hat , trace=trace )
}

