
library(MASS)
library(nets)

N <- 10
T <- 100
P <- 2

A <- array(0,dim=c(N,N,P))
C <- matrix(0,N,N)

A[,,1] <- 0.1*diag(N)
A[,,2] <- 0.5*diag(N)

#
A[1,3,1] <- 0.2
A[2,3,2] <- 0.4

y      <- matrix(0,T,N)
eps    <- mvrnorm(T,rep(0,N),0.1*diag(1,N))

for( t in (P+1):T ){
  for( l in 1:P ){
    y[t,] <- y[t,] + A[,,l] %*% y[t-l,]
  }
  y[t,] <-  y[t,] +  eps[t,]
}

alpha  <- c()
for( p in 1:P ){
  for( i in 1:N ){
    alpha <- c(alpha,A[i,,p])
  }
}

#
lambda  <-1
results <- nets(y,P,lambda=lambda,CN=FALSE,verbose=TRUE) 

if( N < 6 ){
  for( p in 1:P ){ cat('A',p,'\n'); print( round(cbind(A[,,p],rep(NA,N),results$A.hat[,,p]),2) )}
}

results$g.adj

