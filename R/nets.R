
.packageName <- "nets"

.First.lib <- function(lib, pkg){ library.dynam("nets", pkg, lib) }

.onLoad.lib <- function(lib, pkg){ library.dynam("nets", pkg, lib) }

nets <- function( y , type='lrpc' , algorithm='default' , p=1 , lambda=NULL , verbose=FALSE ){ 

	# control for errors
	
	# redirect to the appropriate network estimation routine
	network <- switch( type ,
		lrpc=nets.lrpc( y , p , lambda , verbose ),
		pc=nets.pc( y , lambda , verbose ),
		g=nets.g( y , p , lambda , verbose ) )

	# prepare igraph stuff
	if( type=='lrpc' | type=='pc' ){
		ig <- graph.adjacency( network$Adj , mode='undirected')
	}
	else {
		ig <- graph.adjacency( network$Adj , mode='directed')
	}
	V(ig)$label <- V(ig)$name	
	network$ig <- ig

	# type
	network$type <- type
	network$T    <- nrow(y)
	network$N    <- ncol(y)

	class(network) <- 'nets'

	network
}

plot.nets <- function( network , ... ){ 
  	plot( network$ig , ... )
}

print.nets <- function( network ) {

	labels <- list( lrpc='Long Run Partial Correlation' , pc='Partial Correlation' , g='Granger' )

	if( network$type=='g' ){ 
		nedge <- sum( network$Adj ); 
		medge <- network$N^2 
	}  
	else { 
		nedge <- sum( network$Adj[ lower.tri(network$Adj) ]!=0 ); 
		medge <- ((network$N-1)*network$N)/2.0
	}

	cat( ' ' , labels[[ network$type ]] , ' Network\n' , sep='' )
	cat( ' Number of Detected Edges: ' , nedge , ' (' , round( (nedge/medge)*100 , 1 ) , '%)\n' , sep='' )

}

# Long Run Partial Correlation Network
nets.lrpc <- function(y,lambda,verbose){

	T <- nrow(y)
	N <- ncol(y)
	labels <- dimnames(y)[[2]]

	# G  
	network.G <- nets.g(y,p=p,lambda=lambda[1],verbose=verbose)

	# C
	eps <- results.G$eps 
	results.C <- nets.cov.sfit.alt(eps,lambda=lambda[2],verbose=verbose)
	K <- results.C$K

	# put results together
	G <- matrix(0,N,N)
	G[ col(G)==row(G) ] <- 1
	G <- G - network$G

	KLR <- t(G) %*% K %*% G
	Adj <- (KLR != 0)*1
	Adj[ row(Adj.lr)==col(Adj.lr) ] <- 0

	# packaging results
	dimnames(KLR) <- list(labels,labels)
	dimnames(Adj) <- list(labels,labels)
	dimnames(G)   <- list(labels,labels)
	dimnames(K)   <- list(labels,labels)

}

# Partial Correlation Network
nets.pc <- function(y,lambda,verbose){

	T <- nrow(y)
	N <- ncol(y)
	labels <- dimnames(y)[[2]]
	if( is.null(labels) ) labels <- paste('V',1:N,sep='') 

	results <- nets.space(y,lambda,verbose)	
	bic <- 0

	C <- results$K 
	Adj <-(C!=0)*1
	Adj[ row(Adj)==col(Adj) ] <- 0

	dimnames(C)   <- list(labels,labels)
	dimnames(Adj) <- list(labels,labels)

	network     <- list()
	network$C   <- C
	network$bic <- bic
	network$Adj <- Adj

	network
}

# Granger Network
nets.g <- function(y,p=1,lambda=0,v=NULL,w='adaptive',verbose=FALSE,procedure='shooting'){

	T <- nrow(y)
	N <- ncol(y)
	labels <- dimnames(y)[[2]]
	if( is.null(labels) ) labels <- paste('V',1:N,sep='')
	
	bic <- 0
	A   <- array( 0 , dim=c(p,N,N) )
	G   <- matrix( 0 , N , N )

	for( i in 1:N )
	{
		Y <- matrix( 0 , T-p , 1 ) 
		X <- matrix( 0 , T-p , p*N )
		
		Y[] <- y[ (p+1):T , i]
		for( l in 1:p )
		{
			X[, ( (p-1)*N + p ):( p*N ) ] <- y[ (p+1-l):(T-l) , ]
		}	
	
		results <- nets.alasso( Y , X , w=w, lambda=lambda , verbose=verbose , procedure=procedure )

		for( l in 1:p )
		{
			A[l,i,] <- results$theta[ ( (p-1)*N + p ):( p*N ) ]
		}

		bic <- bic + T * log( sum(eps^2) ) + log(T) * sum( results$theta != 0 ); 
	}
	for( i in 1:p )
	{
		G <- G + A[l,,]
	}

	Adj <-(G!=0)*1

	dimnames(G)   <- list(labels,labels)
	dimnames(Adj) <- list(labels,labels)

	network     <- list()
	network$A   <- A
	network$G   <- G
	network$bic <- bic
	network$Adj <- Adj

	network
}

nest.pc.search <- function(y){}

nest.lrpc.search <- function(y){}

nest.g.search <- function(y){}

# Adaptive Lasso
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

# SPACE Algorithm
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

