
library(nets)

# parameters
N <- 10
T <- 100
P <- 3

A <- array(0,dim=c(N,N,P))

A[3,7,1] <- 0.2
A[4,1,1] <- 0.2
A[9,1,2] <- 0.2
A[2,5,3] <- 0.2
A[8,3,3] <- 0.2

g.adj  <- A[,,1]!=0 | A[,,2]!=0 | A[,,3]!=0

# simulate var
y      <- matrix(0,T,N)
eps    <- matrix(rnorm(T*N),T,N)

for( t in (P+1):T ){
  for( l in 1:P ){
    y[t,] <- y[t,] + A[,,l] %*% y[t-l,]
  }
  y[t,] <-  y[t,] + eps[t,]
}

matplot(y,t='l')

# estimate var
lambda  <- 10
results <- nets(y,P,lambda=lambda,CN=FALSE,verbose=TRUE) 

g.adj.hat <- results$g.adj

# plot the network
granger.network <- graph.adjacency( g.adj.hat , mode='directed' )

plot( granger.network )
