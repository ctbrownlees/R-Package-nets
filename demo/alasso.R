
# The alasso demo simulates a sparse regression model and estimates the coefficient by
# least squares and by the adaptive LASSO for a grid of tuning lasso parameters.

# Parameters
N    <- 100
P    <- 10
P.NZ <- 0.2

# Simulate Sparse Regression Model
nonzero <- rbinom(P,1,P.NZ)
theta.true <- rep( 0 , P );
theta.true[ nonzero==1 ] <- rnorm( sum(nonzero) , 0 , 1)

X <- matrix( rnorm(N*P,0,1) , N , P )
X[,1] <- 1
 
y <- X %*% theta.true + rnorm(N,0,1)

#par( mfrow=c(5,2) )
#for( i in 1:P ){ plot(y,X[,i]) }

readline("\nType  <Return>\t to Continue") 

# LEAST SQUARES
theta.ls <- as.vector(coef(lm.lasso <- lm( y ~ 0+X )))

# LASSO
# notice: range begins with large lambda and then decreases
lambda.range <- rev( c(seq(0,10,1),100) )
L <- length(lambda.range)

theta.lasso   <- matrix( 0 , P , L )
theta.sig2err <- rep( 0 , L )
for( i in 1:L )
{
	results <- alasso( y , X , lambda=lambda.range[i] , procedure="shooting" )

	theta.lasso[,i] <- results$theta
	theta.sig2err[i] <- results$sig2err
}
readline("\nType  <Return>\t to Continue") 

# print
print( cbind( theta.lasso[,c(1,length(lambda.range))] ,  theta.ls , theta.true ) )

#
trace.theta  <- theta.lasso 
trace.lambda <- lambda.range
for( i in 1:L ) trace.theta[,i] <- trace.theta[,i]/max(abs(theta.ls))
trace.theta <- t(trace.theta)

matplot( 1/(1+lambda.range) , trace.theta , t='b' , lwd=3 , ylab='coefficients' , xlab='1/(1+lambda)' , main="Coefficient Trace" )
grid()

readline("\nType  <Return>\t to Continue") 

#
plot( 1/(1+lambda.range) , 1-theta.sig2err/var(y) , t='h' , lwd=4 , col='darkred' , ylab='R Squared' , xlab='1/(1+lambda)' , main="R2 Trace" )
grid()



