
library(MASS)
library(nets)

#set.seed(123456)

N <- 100
T <- 200

C      <- diag(N)
C[3,4] <- -0.4
C[4,3] <- C[3,4]

#
PC    <- -diag( diag(C)**(-0.5) ) %*% C %*% diag( diag(C)**(-0.5) )
Sig    <- solve(C)

y    <- mvrnorm(T,rep(0,N),Sig)

#
lambda  <- 0.5
system.time( results1 <- nets(y,lambda=lambda*T,verbose=TRUE,GN=FALSE,algorithm='activeshooting') )

C.hat.1 <- results1$C.hat

library(space)
results <- space.joint(y,lam1=lambda*T,iter=1)
C.hat.2   <- -diag(sqrt(results$sig.fit)) %*% results$ParCor %*% diag(sqrt(results$sig.fit))
for( i in 1:N ){ C.hat.2[i,i] <- -C.hat.2[i,i] }

print( round( cbind(C,rep(NA,N),C.hat.1,rep(NA,N),C.hat.2) , 2 ) )
