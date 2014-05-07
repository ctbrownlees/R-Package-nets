
require(nets)

library(MASS)

# Problem Parameters
N   <- 10
T   <- 500
P.NZ <- 0.3

# Simulations
S   <- 10

lambda.range <- c(1,2,3,4,5,10,20,100,200,300,500,750,1000,2000)
nlambda      <- length(lambda.range)

mse    <- matrix( 0 , S , nlambda )
tp     <- matrix( 0 , S , nlambda )
tn     <- matrix( 0 , S , nlambda )
fp     <- matrix( 0 , S , nlambda )
fn     <- matrix( 0 , S , nlambda )

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
	
	for( i in 1:nlambda ){

		network <- nets( y, type='pc' , lambda=lambda.range[i] )

		# Some Metrics
		mse[s,i] <- norm( SigInv - network$C , type='F' )
		tp[s,i] <- mean( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]!=0 ] != 0 )
		tn[s,i] <- mean( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]==0 ] == 0 )
		fp[s,i] <- mean( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]==0 ] != 0 )
		fn[s,i] <- mean( (network$C[lower.tri(network$C)])[ SigInv[lower.tri(SigInv)]!=0 ] == 0 )
	}
}
cat('\n')

# More Metrics
tpr <- rep(0,nlambda)
fpr <- rep(0,nlambda)

for( i in 1:nlambda ){
	tpr[i] <- mean( tp[,i] / (tp[,i] + fn[,i]) )
	fpr[i] <- mean( fp[,i] / (fp[,i] + tn[,i]) )
}

plot( fpr , tpr , t='l' , main='ROC' , col='darkblue' , lwd=3 , ylim=c(0,1) , xlim=c(0,1))

plot( lambda.range , colMeans(mse) , t='l' , main='MSE' , col='darkblue' , lwd=3 )

