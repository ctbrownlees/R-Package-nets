
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
	
	# input fix
	# sort & unique  
	lambda <- unique( sort( lambda ) )
	
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
	network$type      <- type
	network$criterion <- select
	network$y         <- y
	network$T         <- nrow(y)
	network$N         <- ncol(y)
	network$L         <- length(lambda)

	class(network) <- 'nets'

	# sanity checks
	if( length(lambda) > 0 )
	{
		# check 1: min lambda shrinks everyone to zero?
		#if( sum(network$theta[, network$lambda==max(network$lambda) ])==0 )
		#{
		#}

		# check 2: max lambda shrinksd does not shrink anyone?

		# check 3: optimal lambda is on the boundary?
	}

	network
}

plot.nets <- function( x , ... ){ 

	# check what to plot
	args  <- list(...)
	if( "what" %in% names(args) ) what <- args['what']
	else what='graph'

	# error check
	if( !any( what==c('graph','ctrace','r2trace') ) ){
		stop("The 'what' parameter has to be set to either 'graph','ctrace' or 'r2trace'")
	}

	if( what=='graph' ){
		# ugh! 
		plot( x$ig )
	}
	if( what=='ctrace' ){
		lambda <- 1/(1+network$lambda[ rank(-network$lambda)]) 
		ctrace <- t(network$theta)[ rank(-network$lambda), ] / max(network$theta[,1])
		lambda.hat <- 1/(1+network$lambda[ network$select ])

		# plot
		matplot(lambda,ctrace,t='b',pch=1,lwd=2,main='C-Trace',xlab='1/(1+lambda)',ylab='standardised coefficients',xpd=FALSE)
		abline(v=lambda.hat,col='gray',lwd=5,xpd=FALSE)	
	}
	if( what=='trace' ){
		# TODO: trace of information criteria
	}
	if( what=='r2trace' ){
		r2trace  <- 100*(1-network$sig2err/mean(diag(var(network$y))))[ rank(-network$lambda) ]

		color    <- rep( 'lightblue' , network$L )
		color[ rev(network$select) ] <- 'gray'

		# plot
		barplot( r2trace , lwd=2 , main='R2-Trace' , xlab='1/(1+lambda)' , ylab='R2' , col=color , yaxt="n" , xaxt="n" )
		par(las=2)
		axis(2, axTicks(2), labels=sprintf("%.1f%%", axTicks(2)) )
		box()
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

	cat( ' Network Type: ' , labels[[ x$type ]] , '\n' , sep='' )
	cat( ' Time Series Dimension: T=',network$T,' N=',network$N,'\n',sep='')
	if( network$L > 1) cat( ' Regularization Method: ',network$criterion,'\n',sep='')
	cat( '\n')
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

	crit    <- matrix(0,1,2,dimnames=list(NULL,c('aic','bic'))) 
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
		eps[,i]       <- results$eps
		crit[1,'bic'] <- crit[1,'bic'] + T*log( (1/T)*sum(eps[,i]^2) ) + sum( results$theta !=0 )*log(T)
		crit[1,'aic'] <- crit[1,'aic'] + T*log( (1/T)*sum(eps[,i]^2) ) + sum( results$theta !=0 )
		sig2err[i]    <- results$sig2err
	}
	cat('\n')

	for( i in 1:p ) {
		G <- G + A[l,,]
	}

	Adj <-(G!=0)*1
	Adj[ row(Adj)==col(Adj) ] <- 0

	#dimnames(A)   <- list(labels,labels)
	dimnames(G)   <- list(labels,labels)
	dimnames(Adj) <- list(labels,labels)

	# assembling the network object
	# common stuff
	network         <- list()
	network$names   <- labels
	network$lambda  <- lambda
	network$select  <- 1 
	network$theta   <- A[1:(N*N*p)] 
	network$sig2err <- mean(sig2err)
	network$Adj     <- Adj
	network$crit    <- crit 

	# granger stuff
	network$granger     <- list()
	network$granger$p   <- p
	network$granger$A   <- A
	network$granger$G   <- G
	network$granger$eps <- eps

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

.nets.g.search <- function(y,p,lambda,select,verbose){

	L <- length(lambda)

	# if there is only one lambda, no search
	if( L == 1 ) return( .nets.g(y,p,lambda,verbose) )

	# otherwise, let's see which is the (max) lambda tha achieves the min 
	trace <- rep(0,L)
	networks <- list()
	
	theta   <- matrix( 0 , ncol(y)*ncol(y)*p , L )
	crit    <- matrix( 0,L,2,dimnames=list(NULL,c('aic','bic')))
	sig2err <- rep( 0 , L ) 

	for( i in 1:L ) {
		networks[[i]] <- .nets.g(y,p,lambda[i],verbose)
		trace[i] <- networks[[i]]$crit[ 1, select ]

		theta[,i] <- networks[[i]]$theta
		crit[i,]  <- networks[[i]]$crit
		sig2err[i] <- networks[[i]]$sig2err
	}
	idx <- max( (1:L)[ min(trace)==trace ] )

	# assembling the network object
	# common stuff
	network         <- list()
	network$lambda  <- lambda
	network$select  <- (1:L)==idx
	network$theta   <- theta
	network$crit    <- crit
	network$sig2err <- sig2err
	network$Adj     <- networks[[ idx ]]$Adj

	# granger stuff 
	network$granger     <- list()
	network$granger$A   <- networks[[ idx ]]$granger$A
	network$granger$G   <- networks[[ idx ]]$granger$G
	network$granger$eps <- networks[[ idx ]]$granger$eps
	
	network
}

# Adaptive Lasso
.nets.alasso <- function(y,X,lambda,maxiter=100,w='adaptive',verbose=FALSE,procedure='activeshooting'){
	
	M <- nrow(y)
	N <- ncol(X)

	toll    <- 1e-6
	
	# check inputs
	if( any( !is.finite(y) ) ){ stop('The response vector contains non finite values.') }
	if( any( !is.finite(X) ) ){ stop('The data matrix contains non finite values.') }
	if( lambda < 0 ){ stop('The ALASSO penalty is negative') }

	# adaptive lasso weights
	if( w=='adaptive' ) {
		if( ncol(X) < nrow(y) ){
			# beta.pre <- coef( lm( y ~ 0+X ) )
			# w <- 1/abs(beta.pre)
			w <- rep(1,N)
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
        results <- .C(procedure, theta=as.double(rep(0,N)) , as.double(y) , as.double(X) , as.double(lambda) , as.double(w) , as.double(theta.init) , as.integer(M) , as.integer(N) , as.integer(verbose) , as.integer(init), as.integer(maxiter) , PACKAGE="nets" )

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

	print( results$ivar )

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

