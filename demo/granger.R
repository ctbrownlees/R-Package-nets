
library(nets)

# parameters
N <- 10
T <- 1000
P <- 3

A <- array(0,dim=c(N,N,P))

A[3,7,1] <- 0.2
A[4,1,1] <- 0.2
A[9,1,2] <- 0.2
A[2,5,3] <- 0.2
A[8,3,3] <- 0.2

g.adj  <- 1 * ( A[,,1]!=0 | A[,,2]!=0 | A[,,3]!=0 )

# simulate var
y      <- matrix(0,T,N)
mu     <- matrix(0,T,N)
eps    <- matrix(rnorm(T*N),T,N)

for( t in (P+1):T ){
  for( l in 1:P ){
    mu[t,] <- mu[t,] + A[,,l] %*% y[t-l,]
  }
  y[t,] <-  mu[t,] + eps[t,]
}

matplot(y,t='l')

# estimate var
lambda  <- 14
mdl <- nets(y[1:round(0.9*T),],P,lambda=lambda,CN=FALSE,verbose=TRUE) 

g.adj.hat <- mdl$g.adj

g.adj.v     <- g.adj[1:(N*N)]
g.adj.hat.v <- g.adj.hat[1:(N*N)] 

tpr <- mean( g.adj.hat.v[ g.adj.v!=0 ] != 0 )
fdr <- mean( g.adj.hat.v[ g.adj.v==0 ] != 0 )

# plot the network
granger.network     <- graph.adjacency( g.adj , mode='directed' )
granger.network.hat <- graph.adjacency( g.adj.hat , mode='directed' )

par( mfrow=c(1,2) )
plot( granger.network )
plot( granger.network.hat )

pred <- predict(mdl,y[round(0.9*T+1):T,])

par( mfrow=c(1,1) )
i <- 3

matplot( cbind( mu[round(0.9*T+1):T,i] , pred$y.hat[,i] ) )
