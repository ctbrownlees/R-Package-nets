
library(MASS)
library(nets)

set.seed(12345)

N <- 10
P <- 2
L <- 40
T <- 1000

A <- array(0,dim=c(N,N,P))
C <- matrix(0,N,N)

A[,,1]   <- 0.7 * diag(N) 
A[,,2]   <- 0.2 * diag(N) 
A[1,2,1] <- 0.2
A[4,3,2] <- 0.2

C      <- diag(N)
C[1,1] <- 2
C[4,2] <- -0.2
C[2,4] <- -0.2
C[1,3] <- -0.1
C[1,3] <- -0.1

g.adj <- 1*( A[,,1]!=0 | A[,,2]!=0 ) - diag(N)
c.adj <- 1*(C!=0) - diag(N)

Sig    <- solve(C)

alpha  <- c()
for( p in 1:P ){
  for( i in 1:N ){
    alpha <- c(alpha,A[i,,p])
  }
}

c      <- diag(C)
PC     <- -diag( c**(-0.5) ) %*% C %*% diag( c**(-0.5) )
rho    <- PC[ upper.tri(PC) ]

y      <- matrix(0,T,N)
eps    <- mvrnorm(T,rep(0,N),Sig)

for( t in (P+1):T ){
  for( l in 1:P ){
    y[t,] <- y[t,] + A[,,l] %*% y[t-l,]
  }
  y[t,] <-  y[t,] +  eps[t,]
}

#
net <- nets(y,P,lambda=2,verbose=TRUE) 

cbind( net$g.adj , g.adj )

cbind( net$c.adj , c.adj )
