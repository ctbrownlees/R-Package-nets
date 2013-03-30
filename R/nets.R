
.packageName <- "nets"

nets <- function( y , type='lrpc' , algorithm='default' , p=1 , lambda=stop("shrinkage parameter 'lambda' has not been set") , select='bic' , verbose=FALSE ){ 

	# input check
	if( !any( type==c('lrpc','pc','g') ) ){
		stop("The 'type' parameter has to be set to either 'lrpc', 'pc' or 'g'")
	}
	if( !is.data.frame(y) & !is.matrix(y) ){
		stop("The 'y' parameter has to be a TxN matrix or a data.frame of data")
	}
	if( p<0 ){
		stop("The 'p' parameter has to be nonnegative")
	}
	if( !any( select==c('aic','bic') ) ){
		stop("The 'type' parameter has to be set to either 'aic' or 'bic'")
	}
	if( !is.logical(verbose) ){
		stop("The 'verbose' parameter has to be TRUE or FALSE")
	}

	# std data (this should be optional)
	for( i in 1:ncol(y) ) y[,i] <- (y[,i]-mean(y[,i]))/sd(y[,i])
	
	# redirect to the appropriate network estimation routine
	network <- switch( type ,
		lrpc=.nets.lrpc.search( y , p , lambda , select , verbose ),
		pc=.nets.pc.search( y , lambda , select , verbose ),
		g=.nets.g.search( y , p , lambda , select , verbose ) )

	# prepare the igraph stuff
	if( type=='lrpc' | type=='pc' ){
		ig <- graph.adjacency( network$Adj , mode='undirected')
	}
	else {
		ig <- graph.adjacency( network$Adj , mode='directed')
	}
	V(ig)$label <- V(ig)$name	
	network$ig <- ig

	# network characteristics 
	network$type <- type
	network$T    <- nrow(y)
	network$N    <- ncol(y)

	class(network) <- 'nets'

	network
}

plot.nets <- function( x , ... ){ 

	if( !exists("what") ) {
		plot( x$ig , ... )
		return()
	}

	if( what=='r2trace' )
	{
	}

	if( what=='ctrace' )
	{
	}

}

print.nets <- function( x , ... ) {

	labels <- list( lrpc='Long Run Partial Correlation' , pc='Partial Correlation' , g='Granger' )

	if( x$type=='g' ){ 
		nedge <- sum( x$Adj ); 
		medge <- x$N^2 
	}  
	else { 
		nedge <- sum( x$Adj[ lower.tri(x$Adj) ]!=0 ); 
		medge <- ((x$N-1)*x$N)/2.0
	}

	cat( ' ' , labels[[ x$type ]] , ' Network\n' , sep='' )
	cat( ' Number of Detected Edges: ' , nedge , ' (' , round( (nedge/medge)*100 , 1 ) , '%)\n' , sep='' )
}

# Long Run Partial Correlation Network
.nets.lrpc <- function(y,p,lambda,verbose){

	y <- as.matrix(y)
	T <- nrow(y)
	N <- ncol(y)
	labels <- dimnames(y)[[2]]

	# G  
	network.G <- .nets.g(y,p=p,lambda=lambda[1],verbose=verbose)

	# C
	network.C <- .nets.pc(network.G$eps,lambda=lambda[2],verbose=verbose)

	# put results together
	C <- network.C$C
	G <- matrix(0,N,N)
	G[ col(G)==row(G) ] <- 1
	G <- G - network.G$G

	KLR <- t(G) %*% C %*% G
	Adj <- (KLR != 0)*1
	Adj[ row(Adj)==col(Adj) ] <- 0

	LRPC  <- diag(1/sqrt(diag(KLR))) %*% KLR %*% diag(1/sqrt(diag(KLR)))
	LRPC[ row(LRPC)!=col(LRPC) ] <- -LRPC[ row(LRPC)!=col(LRPC) ] 

	# packaging results
	dimnames(KLR) <- list(labels,labels)
	dimnames(Adj) <- list(labels,labels)
	dimnames(G)   <- list(labels,labels)
	dimnames(KLR) <- list(labels,labels)

	network      <- list()
	network$C    <- C
	network$G    <- G
	network$KLR  <- KLR
	network$Adj  <- Adj
	network$LRPC <- LRPC 

	network
}

# Partial Correlation Network
.nets.pc <- function(y,lambda,verbose){

	y <- as.matrix(y)
	T <- nrow(y)
	N <- ncol(y)
	labels <- dimnames(y)[[2]]
	if( is.null(labels) ) labels <- paste('V',1:N,sep='') 

	results <- .nets.space(y,lambda,verbose)	
	crit <- list( bic=0 , aic=0 )	

	C   <- results$K 
	Adj <-(C!=0)*1
	Adj[ row(Adj)==col(Adj) ] <- 0
	PC  <- diag(1/sqrt(diag(results$K))) %*% results$K %*% diag(1/sqrt(diag(results$K)))
	PC[ row(PC)!=col(PC) ] <- -PC[ row(PC)!=col(PC) ]

	# measure goodness of fit
	# TODO: INNEFICIENT!!! Check the matrix algebra to simplify this!
	BETA <- matrix(0,N,N)
	for( i in 1:N ){
		for( j in 1:N ){
			BETA[i,j] = -C[i,j]/C[i,i]
		}
	} 
	BETA[ row(BETA)==col(BETA) ] <- 0
	for( i in 1:N ){
		eps <- -y %*% BETA[,i]
		crit$bic <- crit$bic + T*log( (1/T)*sum(eps^2) ) + sum( results$theta !=0 )*log(T)
		crit$aic <- crit$aic + T*log( (1/T)*sum(eps^2) ) + sum( results$theta !=0 )
	}

	dimnames(C)   <- list(labels,labels)
	dimnames(Adj) <- list(labels,labels)
	dimnames(PC)  <- list(labels,labels)

	network      <- list()
	network$C    <- C
	network$crit <- crit
	network$Adj  <- Adj
	network$PC   <- PC 

	network
}

# Granger Network
.nets.g <- function(y,p=1,lambda=0,v=NULL,w='adaptive',verbose=FALSE,procedure='shooting'){

	y <- as.matrix(y)
	T <- nrow(y)
	N <- ncol(y)
	labels <- dimnames(y)[[2]]
	if( is.null(labels) ) labels <- paste('V',1:N,sep='')

	crit    <- list( bic=0 , aic=0 )	
	A       <- array( 0 , dim=c(p,N,N) )
	G       <- matrix( 0 , N , N )
	eps     <- matrix( 0 , T-p , N )
	sig2err <- rep(0,N)

	for( i in 1:N )	{
		# prepare stuff
		Y <- matrix( 0 , T-p , 1 ) 
		X <- matrix( 0 , T-p , p*N )
		
		Y[] <- y[ (p+1):T , i]
		for( l in 1:p ){
			X[, ( (p-1)*N + p ):( p*N ) ] <- y[ (p+1-l):(T-l) , ]
		}
		results <- .nets.alasso( Y , X , w=w, lambda=lambda , verbose=verbose , procedure=procedure )

		# measure goodness of fit
		for( l in 1:p ){
			A[l,i,] <- results$theta[ ( (p-1)*N + p ):( p*N ) ]
		}	
		eps[,i]     <- results$eps
		crit$bic    <- crit$bic + T*log( (1/T)*sum(eps[,i]^2) ) + sum( results$theta !=0 )*log(T)
		crit$aic    <- crit$aic + T*log( (1/T)*sum(eps[,i]^2) ) + sum( results$theta !=0 )
		sig2err[i]  <- results$sig2err
	}
	cat('\n')

	for( i in 1:p ) {
		G <- G + A[l,,]
	}

	Adj <-(G!=0)*1
	Adj[ row(Adj)==col(Adj) ] <- 0

	dimnames(G)   <- list(labels,labels)
	dimnames(Adj) <- list(labels,labels)

	network         <- list()
	network$A       <- A
	network$G       <- G
	network$crit    <- crit
	network$Adj     <- Adj
	network$eps     <- eps

	network$sig2err <- mean(sig2err)
	network$theta   <- A[1:(N*N*p)] 

	network
}

.nets.pc.search <- function(y,lambda,select,verbose){

	if( length(lambda) == 1 ) return( .nets.pc(y,lambda,verbose) )

	trace <- rep(0,length(lambda))
	networks <- list()

	for( i in 1:length(lambda) ) {
		networks[[i]] <- .nets.pc(y,lambda[i],verbose)
		trace[i] <- networks[[i]]$crit[[ select ]]
	}
	idx <- max( (1:length(lambda))[ min(trace)==trace ] )

	# package results
	network <- networks[[ idx ]]
	network$idx   <- idx
	network$trace <- trace

	network
}

.nets.lrpc.search <- function(y,p,lambda,select,verbose){

}

.nets.g.search <- function(y,p,lambda,crit,select,verbose){

	# if there is only one lambda, no search
	if( length(lambda) == 1 ) return( .nets.g(y,p,lambda,verbose) )

	# otherwise, let's see which is the (max) lambda tha achieves the min 
	trace <- rep(0,length(lambda))
	networks <- list()

	theta   <- matrix( 0 , ncol(y)*ncol(y)*p , length(lambda) )
	sig2err <- rep( 0 , length(lambda) ) 

	for( i in 1:length(lambda) ) {
		networks[[i]] <- .nets.g(y,p,lambda[i],verbose)
		trace[i] <- networks[[i]]$crit[[ select ]]

		theta[,i] <- networks[[i]]$theta
		sig2err[i] <- networks[[i]]$sig2err
	}
	idx <- max( (1:length(lambda))[ min(trace)==trace ] )

	# package results
	network         <- networks[[ idx ]]
	network$select  <- idx
	network$trace   <- trace
	network$theta   <- theta
	network$sig2err <- sig2err

	network
}

# Adaptive Lasso
.nets.alasso <- function(y,X,lambda,w='adaptive',verbose=FALSE,procedure='activeshooting'){
	
	M <- nrow(y)
	N <- ncol(X)

	toll    <- 1e-6
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
        results <- .C(procedure, theta=as.double(rep(0,N)) , as.double(y) , as.double(X) , as.double(lambda) , as.double(w) , as.double(theta.init) , as.integer(M) , as.integer(N) , as.integer(verbose) , as.integer(init) , PACKAGE="nets" )

	# packaging results
	results <- list( theta=results$theta , eps=(y-X%*%results$theta) , sig2err=mean((y-X%*%results$theta)**2) )
}

alasso <- .nets.alasso

# SPACE Algorithm
.nets.space <- function(y,lambda,verbose=FALSE)
{
	M <- nrow(y)
	N <- ncol(y)
	
	results <- .C('space', theta=as.double(rep(0,N*(N-1)/2.0)) , ivar=as.double(rep(0,N)) , as.double(y), as.double(lambda), as.integer(M) , as.integer(N) , as.integer(verbose) , PACKAGE="nets" )

	# packaging results
	K <- matrix( 0 , N , N )
	diag(K) <- results$ivar
	for( i in 2:N ){
		for( j in 1:(i-1) ){
			K[i,j] <- -results$theta[ (i-1)*(i-2)/2+j  ]*results$ivar[i]
			K[j,i] <- K[i,j]
		}
	}

	list( K=K , results=results )
}

