
require(nets)

# Problem Parameters
n   <- 1000
p   <- 20
pnz <- 0.2
sigerr <- 1;

# Simulate Regression Model
nonzero <- rbinom(p,1,pnz)
theta <- rep( 0 , p );
theta[ nonzero==1 ] <- rnorm( sum(nonzero) , 0 , 1)

X <- matrix( rnorm(n*p,0,1) , n , p )
X[,1] <- 1 
y <- X %*% theta + rnorm(n,0,sigerr)

results <- nets.alasso( y , X , 4 )

cbind( theta , results$theta )

