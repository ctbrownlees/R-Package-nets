
# Clean up
rm( list=ls() )

# Parameters
S    <- 10
N    <- 100
P    <- 70
P.NZ <- 0.2

# Simulate Sparse Regression Model
nonzero <- rbinom(P,1,P.NZ)
theta.true <- rep( 0 , P );
theta.true[ nonzero==1 ] <- 1 

# Generate Regressors
X <- matrix( rnorm(N*P,0,1) , N , P )
X[,1] <- 1

# Lambda Range
lambda.range <- rev( c( 0 , seq(0.1,10,0.1) , 100 ) )
L <- length(lambda.range)

#
mse <- matrix( 0 , S , L )
tp  <- matrix( 0 , S , L )
tn  <- matrix( 0 , S , L )
fp  <- matrix( 0 , S , L )
fn  <- matrix( 0 , S , L )

# Simulation
for( s in 1:S )
{
	y <- X %*% theta.true + rnorm(N,0,1)

	for( l in 1:L )
	{
		cat('.') 

		results <- alasso( y , X , lambda=lambda.range[l] , procedure="shooting")
		theta.lasso <- results$theta

		mse[s,l] <- mean( ( theta.lasso - theta.true )**2  )
		tp[s,l]  <- sum( theta.lasso[ theta.true != 0  ] != 0 )	
		tn[s,l]  <- sum( theta.lasso[ theta.true == 0  ] == 0 )	
		fp[s,l]  <- sum( theta.lasso[ theta.true == 0  ] != 0 )	
		fn[s,l]  <- sum( theta.lasso[ theta.true != 0  ] == 0 )	
	}
}

# compuet MSE
mse <- colSums(mse)

plot( lambda.range , mse , col='darkblue' , lwd=3)

# Compute and plot the ROC 
tp <- colSums(tp)
tn <- colSums(tn)
fp <- colSums(fp)
fn <- colSums(fn)

tpr <- tp / (tp + fn) 
fpr <- fp / (fp + tn)

plot( fpr , tpr , t='l' , main='ROC' , col='darkblue' , lwd=3 , ylim=c(0,1) , xlim=c(0,1))
grid()


