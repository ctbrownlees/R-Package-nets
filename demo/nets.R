
library(MASS)
library(nets)

set.seed(12345)

N <- 5
P <- 2
L <- 40
T <- 200

A <- array(0,dim=c(N,N,P))
C <- matrix(0,N,N)

A[,,1]   <- 0.7 * diag(N) 
A[,,2]   <- 0.2 * diag(N) 
A[1,1,1] <- 0

C      <- diag(N)
C[1,1] <- 2
C[4,2] <- -0.2
C[2,4] <- -0.2
C[1,3] <- -0.1
C[1,3] <- -0.1

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
lambda  <- 0.01
results2 <- nets(y,P,lambda=lambda*T,verbose=TRUE) 
