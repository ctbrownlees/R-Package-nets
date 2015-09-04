
.packageName <- "nets"

nets <- function( y , p=1 , GN=TRUE , CN=TRUE , lambda=stop("shrinkage parameter 'lambda' has not been set") , verbose=FALSE ){

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
  if( GN==FALSE & CN==FALSE ){
    stop("At least A or G have to be true")  
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
	if( T>N*P ){
	  if( GN == TRUE ){
	    x.aux <- matrix( 0 , T , N*P )
	    for( p in 1:P ){
	      x.aux[(p+1):T, ((p-1)*N+1):(p*N) ] <- y[1:(T-p),]
	    }
	    
	    reg <- lm( y ~ 0+x.aux )
      A     <- coef(reg)
	    eps   <- reg$resid
	    
	    alpha <- c()
	    for( p in 1:P ){
	      alpha <- c( alpha , ( A[((p-1)*N+1):(p*N),] ) [1:(N*N)] )
	    }
	    alpha.weights <- 1/abs(alpha)
      
	    # warmup
      for( i in 1:N ){
        for( j in 1:N ){
          for( p in 1:P ){
            trim <- cor.test( y[(p+1):T,i] , y[1:(T-p),j] )$p.value < 0.10
            A[(p-1)*N+j,i] <- A[(p-1)*N+j,i] * trim
          }
        }        
      }
	    alpha <- c()
	    for( p in 1:P ){
	      alpha <- c( alpha , ( A[((p-1)*N+1):(p*N),] ) [1:(N*N)] )
	    }
	    
	  }
	  else{
	    eps           <- y
	    alpha         <- c()
	    alpha.weights <- c()
	  }
	  
	  if( CN == TRUE ){
	    C.hat       <- solve( cov(eps) )
	    c.hat       <- diag(C.hat)
	    PC          <- -diag( c.hat**(-0.5) ) %*% C.hat %*% diag( c.hat**(-0.5) )

      rho         <- PC[ lower.tri(PC) ]
      rho.weights <- 1/abs(rho)

      # warmup
      for( i in 2:N ){
        for( j in 1:(i-1) ){
          PC[i,j] <- PC[i,j] * (cor.test(eps[,i],eps[,j])$p.value < 0.1)
        }
      }
	    rho         <- PC[ lower.tri(PC) ]
	  }
	  else{
	    rho         <- c()
	    rho.weights <- c()
	  }
	}
	else{
	  stop("Not implemented yet")
	}	
  
  #print( alpha )
  #print( rho )
	
  # call nets
	run <- .C("nets3",
	          alpha        =as.double(alpha),
	          rho          =as.double(rho), 
	          alpha.weights=as.double(alpha.weights),
	          rho.weights  =as.double(rho.weights),
	          lambda       =as.double(lambda),
	          y            =as.double(y),
	          T            =as.integer(T),
	          N            =as.integer(N),
	          P            =as.integer(P),
	          c.hat        =as.double(c.hat),
	          GN           =as.integer(GN),
	          CN           =as.integer(CN),
	          v            =as.integer(verbose) )
  
  # package results
	obj <- list()
	class(obj)    <- 'nets'
	obj$T         <- T
	obj$N         <- N
	obj$P         <- P
	obj$lambda    <- lambda
	
  if( GN == TRUE ){
  	A.hat <- array(0,dim=c(N,N,P))
  	for( p in 1:P ){
  	  for( i in 1:N ){
  	    A.hat[i,,p] <- run$alpha[ ((p-1)*N*N+(i-1)*N+1):((p-1)*N*N+i*N) ]
	    }
	  }
  	dimnames(A.hat)[[1]] <- labels
  	dimnames(A.hat)[[2]] <- labels    
  	obj$A.hat     <- A.hat 
  	obj$alpha.hat <- run$alpha 
  }
  if( CN == TRUE ){
	  C.hat <- matrix(0,N,N)
    for( i in 1:N ){
      C.hat[i,i] <- run$c[i]
      for( j in setdiff(1:N,i) ){
        c_ij       <- -run$rho[ (max(i,j)-1)*(max(i,j)-2)/2 + min(i,j) ] * sqrt( run$c.hat[i] * run$c.hat[j] )
        C.hat[i,j] <- c_ij
        C.hat[j,i] <- c_ij
      }
    }
	  dimnames(C.hat)[[1]] <- labels
	  dimnames(C.hat)[[2]] <- labels
	  obj$C.hat     <- C.hat 
	  obj$rho.hat   <- run$rho
  }
  
  # networks
  I.n <- diag(N)
  
  if( CN == TRUE ){
    PCN  <- -diag( diag(C.hat)**(-0.5) ) %*% C.hat %*% diag( diag(C.hat)**(-0.5) )
	  PCN[row(I.n) == col(I.n) ] <- 1
    obj$pc.net          <- PCN
  }
  
  if( GN == TRUE ){
    DGN <- matrix(0,N,N)
    for( p in 1:P ){
      DGN <- DGN + A.hat[,,p]
    }
	  DGN[row(I.n) == col(I.n) ] <- 0
    obj$dgranger.net    <- DGN
  }
  
  if( CN == TRUE && GN==TRUE ){
    KL    <- t(I.n-DGN ) %*% C.hat %*% ( I.n-DGN )
    LRPCN <- -diag( diag(KL)**(-0.5) ) %*% KL %*% diag( diag(KL)**(-0.5) )
	  LRPCN[row(I.n) == col(I.n) ] <- 1
	  obj$lrpc.net        <- LRPCN
  }

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
   	cat( ' Time Series Dimension: T=',x$T,' N=',x$N,'\n',sep='')
   	cat( ' VAR Lags P=',x$P,'\n',sep='')
}

plot.nets <- function( x ){
  print('plot')
}

summary.nets <- function( x , ... ) {
  print('summary.nets')
  class(x) <- 'summary.nets'
  return(x)
}

print.summary.nets <- function( x , ... ){
  print('summary nets')
}
