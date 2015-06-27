
.packageName <- "nets"

nets <- function( y , p=1 , lambda=stop("shrinkage parameter 'lambda' has not been set") , verbose=FALSE ){

	# input check
	if( !is.data.frame(y) & !is.matrix(y) ){
		stop("The 'y' parameter has to be a TxN matrix or a data.frame of data")
	}
	if( p<0 ){
		stop("The 'p' parameter has to be nonnegative")
	}
	if( !is.logical(verbose) ){
		stop("The 'verbose' parameter has to be TRUE or FALSE")
	}

  # define variables
	T <- nrow(y)
	N <- ncol(y)
	P <- p
	
  if( !is.null(dimnames(y)[[2]]) ){
    labels <- dimnames(y)[[2]]
  } else {
    labels <- paste('V',1:N,sep='')
  }
	
  # pre-estimation
	x.aux <- matrix( 0 , T , N*P )
	for( p in 1:P ){
	  x.aux[(p+1):T, ((p-1)*N+1):(p*N) ] <- y[1:(T-p),]
	}
	
	reg <- lm( y ~ 0+x.aux )
	
	alpha <- c()
	for( p in 1:P ){
	  alpha <- c( alpha , ( coef(reg)[((p-1)*N+1):(p*N),] ) [1:(N*N)] )
	}
	
	K.hat <- solve( cov(reg$resid) )
	kk    <- diag(K.hat)
	PC    <- -diag( kk**(-0.5) ) %*% K.hat %*% diag( kk**(-0.5) )
	rho   <- PC[ lower.tri(PC) ]
	
	alp.weights <- 1/abs(alpha)
	rho.weights <- 1/abs(rho)
	
  # call nets
	run <- .C("nets",
	          alpha        =as.double(alpha),
	          rho          =as.double(rho), 
	          alpha.weights=as.double(alp.weights),
	          rho.weights  =as.double(rho.weights),
	          lambda       =as.double(lambda),
	          y            =as.double(y),
	          T            =as.integer(T),
	          N            =as.integer(N),
	          P            =as.integer(P),
	          kk           =as.double(kk),
	          v            =as.integer(verbose) )
	
  # package results
	A.hat <- array(0,dim=c(N,N,P))
	for( p in 1:P ){
	  for( i in 1:N ){
	    A.hat[i,,p] <- run$alpha[ ((p-1)*N*N+(i-1)*N+1):((p-1)*N*N+i*N) ]
	  }
	}
	C.hat <- matrix(0,N,N)
  P.hat <- matrix(0,N,N)
  for( i in 1:N ){
   C.hat[i,i] <- run$kk[i]
   for( j in setdiff(1:N,i) ){
       c_ij       <-  -run$rho[ (max(i,j)-1)*(max(i,j)-2)/2 + min(i,j) ] * sqrt( run$kk[i] * run$kk[j] )
       C.hat[i,j] <- c_ij
       C.hat[j,i] <- c_ij
    }
  }
	dimnames(A.hat)[[1]] <- labels
	dimnames(A.hat)[[2]] <- labels
	dimnames(C.hat)[[1]] <- labels
	dimnames(C.hat)[[2]] <- labels
	  
	obj <- list( )
	class(obj)    <- 'nets'
	obj$A.hat     <- A.hat 
  obj$C.hat     <- C.hat 
  obj$alpha.hat <- run$alpha 
  obj$rho.hat   <- run$rho
	
  #
  I.n <- diag(N)
  CPCN  <- -diag( diag(C.hat)**(-0.5) ) %*% C.hat %*% diag( diag(C.hat)**(-0.5) )
	CPCN[row(I.n) == col(I.n) ] <- 1
  DGN <- matrix(0,N,N)
  for( p in 1:P ){
    DGN <- DGN + A.hat[,,p]
  }
	DGN[row(I.n) == col(I.n) ] <- 0
  KL    <- t(I.n-DGN ) %*% C.hat %*% ( I.n-DGN )
  LRPCN <- -diag( diag(KL)**(-0.5) ) %*% KL %*% diag( diag(KL)**(-0.5) )
	LRPCN[row(I.n) == col(I.n) ] <- 1
  
  obj$cpc.net         <- CPCN
  obj$dgranger.net    <- DGN
	obj$lrpc.net        <- LRPCN

	return(obj)
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
