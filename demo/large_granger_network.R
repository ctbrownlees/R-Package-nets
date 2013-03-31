
# Clean Up
rm( list=ls() )

# Libraries
require(nets)

# Parameters
N   <- 25
T   <- 500
p   <- 0.1
a   <- 0.3
 
# Simulate a large dimensional VAR(1)
ER <- erdos.renyi.game(N,p,directed=TRUE)
A  <- a * as.matrix( get.adjacency(ER) )

y <- matrix(0,T,N)
eps <- matrix( rnorm(T*N,0,1) , T , N )
for( t in 2:T ){ y[t,] = A %*% y[t-1,] + eps[t,] }

readline("Type <Return> to Continue") 

# Estimate Network 
cat('Estimating Granger Network...')

network <- nets( y, p=1, type='g', lambda=c( seq(1,20,0.5) ) , std=FALSE )

cat('done!\n')

readline("Type <Return> to Continue") 

