
.packageName <- "nets"

.First.lib <- function(lib, pkg)
{
     library.dynam("nets", pkg, lib)

}

nets <- function(){
     a <- c(1,2,3)
     results <- .C("adding", as.double(a), as.integer(length(a)), ab = double(length(a)))
     print( results )
}

plot.nets <- function(network){

}

nets.alasso <- function(y,X,lambda,w='adaptive'){
	
	toll <- 1e-6
	maxiter <- 20
	T <- length(y)
	P <- ncol(X)

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
			# ridge
			w <- rep(1,P)
		}
	}
	else {
		w <- rep(1,P);
	}

	# call shooting algorithm

	results <- list( theta=theta , eps=(y-X%*%theta) )
}


nets.test <- function()
{
	y <- c(1,2,2)

	X <- matrix( 1:9 , 3 ,3 )

	theta <- c(0,0,0)

        results <- .C("shooting", as.double(v), as.double(A), as.double(r) )

	print( results )

	print( A %*% v )
}
