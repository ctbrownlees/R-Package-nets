library(nets)

# Parameters
N    <- 31
P    <- 30
P.NZ <- 0.2
NbIter <- 100

# select Lambda to test
lambda.range <- rev( seq(0,20,0.2) )
L <- length(lambda.range)



nonzero <- rbinom(P,1,P.NZ)
theta.true <- rep( 0 , P );
theta.true[ nonzero==1 ] <- rnorm( sum(nonzero) , 0 , 1)

#For every simulation we record the number of zero-parameter there is in theta.true and in theta.lambda
NbZeros = matrix(0,NbIter,length(lambda.range))
NbZerostrue = matrix(0,NbIter)


for( i in 1:NbIter )
	{#Initialize theta, X and y
	NbZerostrue[i] <- length(which(theta.true == 0))
	X <- matrix( rnorm(N*P,0,1) , N , P )
	X[,1] <- 1
	y <- X %*% theta.true + rnorm(N,0,1)

	#Lasso algo for every Lambda
	theta.lasso   <- matrix( 0 , P , L )
	theta.sig2err <- rep( 0 , L )
	
	for( j in 1:L )
		{
			results <- alasso( y , X , lambda=lambda.range[j] , procedure="shooting" )

			theta.lasso[,j] <- results$theta
			theta.sig2err[j] <- results$sig2err
			NbZeros[i,j] <- length(which(theta.lasso[,j] == 0))
		}


	}

#Compare true number of zeros and lasso
CompareZeros <- NbZeros
for (i in 1:L)
	{CompareZeros[,i] <- NbZeros[,i]-NbZerostrue}

#Mean discrepencies of nonzero parameters depeding on lambda
Mdiscrep <- abs(colSums (CompareZeros, na.rm = FALSE, dims = 1))/100
plot(lambda.range,Mdiscrep)

