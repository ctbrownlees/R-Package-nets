
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
SigInv <- matrix( 0 , N , N )
for (i in 2:N)
{
	for (j in 1:(i-1))
	{
		SigInv[i,j] <- rbinom(1,1,P.NZ) 
		SigInv[j,i] <- SigInv[i,j]
	}
}
SigInv <- SigInv + diag( colSums(SigInv) + 1)

for (s in 1:S)
{ 
	cat('.')
	y <- mvrnorm(T, rep(0,N) , solve(SigInv) )

	network <- nets( y, type='pc' , lambda=seq(10,90,10) )

	# Some Metrics
	tp[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]!=0 ] != 0 )
	tn[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]==0 ] == 0 )
	fp[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]==0 ] != 0 )
	fn[s] <- sum( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]!=0 ] == 0 )
}
cat('\n')

# More Metrics
tpr <- tp / (tp + fn) 
fpr <- fp / (fp + tn)

plot( fpr , tpr , t='l' , main='ROC' , col='darkblue' , lwd=3 , ylim=c(0,1) , xlim=c(0,1))


