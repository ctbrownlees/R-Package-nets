
.packageName <- "nets"

.First.lib <- function(lib, pkg){ library.dynam("nets", pkg, lib) }

.onLoad.lib <- function(lib, pkg){ library.dynam("nets", pkg, lib) }

nets <- function( y , p=1 , lambda.G=NULL , lambda.G.range=NULL , lambda.C=NULL , lambda.C.range=NULL , verbose=FALSE ){ 

	T <- nrow(y)
	N <- ncol(y)
	labels <- dimnames(y)[[2]]
	
	# if labels are nul... do somethin'!
	# labels <- paste('V',1:N,sep='')

	# G  
	results.G <- nets.var.sfit(y,p=p,lambda=lambda.G,lambda.range=lambda.G.range,verbose=verbose)

	G <- matrix(0,N,N)
	G[ col(G)==row(G) ] <- 1
	for( i in 1:p ){ G <- G - results.G$A[i,,] }

	# C
	eps <- results.G$eps 
	results.C <- nets.cov.sfit.alt(eps,lambda=lambda.C,lambda.range=lambda.C.range,verbose=verbose)
	K <- results.C$K

	# all togher
	K.lr <- t(G) %*% K %*% G
	Adj.lr <- (K.lr != 0)*1
	Adj.lr[ row(Adj.lr)==col(Adj.lr) ] <- 0

	# packaging results
	dimnames(K.lr) <- list(labels,labels)
	dimnames(Adj.lr) <- list(labels,labels)
	dimnames(G) <- list(labels,labels)
	dimnames(K) <- list(labels,labels)

	ig <- graph.adjacency( Adj.lr , mode='undirected')
	V(ig)$label <- V(ig)$name
	
	results = list( K.lr=K.lr , Adj.lr=Adj.lr , G=G , K=K , lambda.G=results.G$lambda , lambda.C=results.C$lambda , ig=ig )
	class(results) <- 'nets'

	results
}

plot.nets <- function( network , ... ){ 
  	plot( network$ig , ... )
}

print.nets <- function( network ) {

	cat( ' Number of Detected Edges: ' , sum( network$K.lr[ lower.tri( network$K.lr ) ]!=0 ) , ' (out of ' , length( network$K.lr[ lower.tri( network$K.lr ) ]!=0 ) , ' possible edges)\n' , sep='' )
	cat( ' Granger Links: ' , sum( network$G[ col(network$G)!=row(network$G) ]!=0 ) , ' Contemporaneous Links: ' , sum( network$K[ lower.tri( network$K ) ]!=0 ) , '\n' , sep='')

}

nets.alasso <- function(y,X,lambda,w='adaptive',verbose=FALSE,procedure='shooting'){
	
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

	theta.init <- rep(0,N)
	init <- 0

	# call shooting algorithm
        results <- .C(procedure, theta=as.double(rep(0,N)) , as.double(y) , as.double(X) , as.double(lambda) , as.double(w) , as.double(theta.init) , as.integer(M) , as.integer(N) , as.integer(verbose) , as.integer(init) )

	# packaging results
	results <- list( theta=results$theta , eps=(y-X%*%results$theta) )
}

nets.space <- function(y,lambda,verbose=FALSE)
{
	M <- nrow(y)
	N <- ncol(y)
	
	results <- .C('space', theta=as.double(rep(0,N*(N-1)/2.0)) , ivar=as.double(rep(0,N)) , as.double(y), as.double(lambda), as.integer(M) , as.integer(N) , as.integer(verbose) )

	K <- matrix( 0 , N , N )
	diag(K) <- results$ivar
	for( i in 2:N ){
		for( j in 1:(i-1) ){
			K[i,j] <- -results$theta[ (i-1)*(i-2)/2+j  ]*results$ivar[i]
			K[j,i] <- K[i,j]
		}
	}

	# packaging resuls
	results <- list( K=K , results=results )
}


# SPARSE VAR ESTIMATION
nets.var.sfit <- function(y,p=1,lambda=0,lambda.range=NULL,v=NULL,w='adaptive',verbose=FALSE,procedure='shooting'){

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
	if(verbose) cat('Sparse VAR estimation: ')
	stuff = list()
	trace = matrix(0,length(lambda.range),N*N*p)
	for( li in 1:length(lambda.range) )
	{
		if(verbose) cat('.')

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
		
			results <- nets.alasso( Y , X , w=w, lambda=lambda.range[li] , verbose=verbose , procedure=procedure )

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
	if(verbose) cat('\n')

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

# SPARSE INVERSE COVARIANCE ESTIMATION Barigozzi Brownlees
nets.cov.sfit <- function(y,lambda=NULL,lambda.range=NULL,w='adaptive',verbose=FALSE){

	# 
	T <- nrow(y)
	N <- ncol(y)

	# lambda range
	if( is.null(lambda.range) )
	{
		lambda.range <- lambda
		ic <- rep(0,1)
	}
	else
	{
		lambda.range <- rev( lambda.range )
		ic <- rep( 0 , length(lambda.range) )
	}

	# ESTIMATE
	if(verbose) cat('Sparse Inverse Covariance estimation: ')
	stuff = list()
	for( li in 1:length(lambda.range) )
	{
		if(verbose) cat('.')

		# updating theta
		results <- nets.space( y , lambda=lambda.range[li] , verbose=verbose )

		stuff[[li]]   <- list()
		stuff[[li]]$K <- results$K
	}
	if(verbose) cat('\n')

	opt <- (1:length(lambda.range))[ rank(ic,ties.method="first")==1 ]
	lambda.hat <- lambda.range[ opt ]
	K <- stuff[[opt]]$K

	# result
	list( K=K, lambda=lambda.hat )
}

nets.cov.sfit.alt <- function(y,lambda=NULL,lambda.range=NULL,w='adaptive',verbose=FALSE,procedure='activeshooting'){

	T <- nrow(y)
	N <- ncol(y)

	# lambda range
	if( is.null(lambda.range) )
	{
		lambda.range <- lambda
		bic <- rep(0,1)
		bic <- rep( 0 , 1 )
		rss <- rep( 0 , 1 )
		nze <- rep( 0 , 1 )
	}
	else
	{
		lambda.range <- rev( lambda.range )
		bic <- rep( 0 , length(lambda.range) )
		rss <- rep( 0 , length(lambda.range) )
		nze <- rep(0,length(lambda.range))
	}

	# ESTIMATE
	if(verbose) cat('Sparse Inverse Covariance estimation: ')
	stuff = list()
	for( li in 1:length(lambda.range) )
	{
		if(verbose) cat('.')

		results   <- space.joint(y,lam1=lambda.range[li])

		BETA <- results$ParCor
		for( i in 1:N ){
			for( j in 1:N ){
				BETA[i,j] <- results$ParCor[i,j] * ( results$sig.fit[j]/results$sig.fit[i] )
			}
		}

		results$K <- diag( sqrt(results$sig.fit) ) %*% results$ParCor %*% diag( sqrt(results$sig.fit) )
		results$K[ col(results$K) != row(results$K) ] <- -results$K[ col(results$K) != row(results$K) ]
	
		bic[li] <- 0
		nze[li] <- 0
		rss[li] <- 0
		for( i in 1:N ){ 
			nzero <- 0
			for( j in (1:N)[(1:N)!=i] ){
				eps <- y[,i]
				for( j in (1:N)[(1:N)!=i] ){
					eps <- eps - BETA[i,j] * y[,j]
					if( BETA[i,j]!=0 ) nzero <- nzero + 1;
				}
				if( results$K[i,j]!=0 )	nzero <- nzero + 1;
			}
			bic[li] <- bic[li] + T * log( sum(eps^2) ) + log(T) * nzero
			rss[li] <- rss[li] + T * log( sum(eps^2) )
			nze[li] <- nze[li] + nzero;
		}
		stuff[[li]] <- list()
		stuff[[li]]$K <- results$K
	}
	if(verbose) cat('\n')

	print( bic )
	print( rss )
	print( nze )
	print( bic + log(T)*nze )

	opt <- (1:length(lambda.range))[ rank(bic,ties.method="first")==1 ]
	print( opt )
	lambda.hat <- lambda.range[ opt ]

	results   <- space.joint(y,lam1=lambda.hat)
	results$K <- diag( sqrt(results$sig.fit) ) %*% results$ParCor %*% diag( sqrt(results$sig.fit) )
	results$K[ col(results$K) != row(results$K) ] <- -results$K[ col(results$K) != row(results$K) ]

	K <- results$K

	# result
	list( K=K, lambda=lambda.hat )
}
