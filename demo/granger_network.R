
library(nets)

N <- 10
T <- 100
P <- 3

A <- array(0,dim=c(N,N,P))

A[,,1] <- 0.1*diag(N)
A[,,2] <- 0.5*diag(N)

y      <- matrix(0,T,N)
eps    <- matrix(rnorm(T*N),T,N)

for( t in (P+1):T ){
  for( l in 1:P ){
    y[t,] <- y[t,] + A[,,l] %*% y[t-l,]
  }
  y[t,] <-  y[t,] + eps[t,]
}

matplot(y,t='l')

#
lambda  <-1
results <- nets(y,P,lambda=lambda,CN=FALSE,verbose=TRUE) 

if( N < 6 ){
  for( p in 1:P ){ cat('A',p,'\n'); print( round(cbind(A[,,p],rep(NA,N),results$A.hat[,,p]),2) )}
}

results$g.adj
