
# Parameters
N    <- 20
P    <- 25
P.NZ <- 0.2
lambda.range <- rev( seq(1,11,5) )
L <- length(lambda.range)

nonzero <- rbinom(P,1,P.NZ)
theta.true <- rep( 0 , P );
theta.true[ nonzero==1 ] <- rnorm( sum(nonzero) , 0 , 1)

	X <- matrix( rnorm(N*P,0,1) , N , P )
	X[,1] <- 1
	y <- X %*% theta.true + rnorm(N,0,1)

	#Lasso algo for every Lambda
	theta.lasso   <- matrix( 0 , P , L )
	theta.sig2err <- rep( 0 , L )

	for( j in 1:L )
	{#plot(lambda.range[j])
	results <- alasso( y , X , lambda=lambda.range[j], 0 , w="non" , procedure="shooting" )
	theta.lasso[,j] <- results$theta
	theta.sig2err[j] <- results$sig2err
	
	}
theta.ls <- as.vector(coef(lm.lasso <- lm( y ~ 0+X )))
cbind(theta.lasso,theta.ls,theta.true)

t(X)%*%y*-2

