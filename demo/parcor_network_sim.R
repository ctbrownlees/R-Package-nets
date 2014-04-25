
require(nets)

library(MASS)

# Problem Parameters
N   <- 10
T   <- 500
P.NZ <- 0.3

# Simulations
S   <- 100

tp     <- matrix( 0 , S , 1 )
tn     <- matrix( 0 , S , 1 )
fp     <- matrix( 0 , S , 1 )
fn     <- matrix( 0 , S , 1 )


# SigInv
nonzero <- rbinom(N*(N-1)/2,1,P.NZ)*runif(N*(N-1)/2)
SigInv <- diag( 1 , N , N )
count <- 1
for (i in 1:(N-1))
{
	for (j in (i+1):N)
	{
		SigInv[i,j] <-  nonzero[count]
		SigInv[j,i] <-  nonzero[count]
		count = count+1	
	}
}


for (s in 1:S)
{ 
	y <- mvrnorm(T, rep(0,N) , solve(SigInv) )

	network <- nets( y, type='pc' , lambda=seq(10,90,10) )



	# Some Metrics
	tp[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]!=0 ] != 0 )
	tn[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]==0 ] == 0 )
	fp[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]==0 ] != 0 )
	fn[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]!=0 ] == 0 )
}

# More Metrics
tpr <- tp / (tp + fn) 
fpr <- fp / (fp + tn)

plot( fpr , tpr , t='l' , main='ROC' , col='darkblue' , lwd=3 , ylim=c(0,1) , xlim=c(0,1))


