
#
require(nets)

# libraries
library(MASS)

mnorm <- function(M){ sqrt( max( eigen( M %*% M )$values ) ) }

# System Dimension
N <- 6
names <- c('A','B','C','D','E','F')

# 
I      <- matrix( 0 , N , N )
A      <- matrix( 0 , N , N , dimnames=list(names,names) )
SigInv <- matrix( 0 , N , N , dimnames=list(names,names) )

# I
I[ row(I)==col(I) ] <- 1

# A
A[1,1] <- 0.71; A[1,2] <- 0.00; A[1,3] <- 0.00; A[1,4] <- 0.00; A[1,5] <- 0.00; A[1,6] <- 0.15;
A[2,1] <- 0.00; A[2,2] <- 0.63; A[2,3] <- 0.00; A[2,4] <- 0.30; A[2,5] <- 0.00; A[2,6] <- 0.00;
A[3,1] <- 0.00; A[3,2] <- 0.00; A[3,3] <- 0.10; A[3,4] <- 0.00; A[3,5] <- 0.00; A[3,6] <- 0.00;
A[4,1] <- 0.00; A[4,2] <- 0.00; A[4,3] <- 0.10; A[4,4] <- 0.22; A[4,5] <- 0.25; A[4,6] <- 0.00;
A[5,1] <- 0.10; A[5,2] <- 0.00; A[5,3] <- 0.00; A[5,4] <- 0.00; A[5,5] <- 0.34; A[5,6] <- 0.00;
A[6,1] <- 0.00; A[6,2] <- 0.00; A[6,3] <- 0.00; A[6,4] <- 0.00; A[6,5] <- 0.00; A[6,6] <- 0.42;

# SigInv
SigInv[1,1] <-  1.00; SigInv[1,2] <-  0.00; SigInv[1,3] <-  0.00; SigInv[1,4] <-  0.00; SigInv[1,5] <- 0.00; SigInv[1,6] <- 0.00;
SigInv[2,1] <-  0.00; SigInv[2,2] <-  1.00; SigInv[2,3] <-  0.00; SigInv[2,4] <-  0.00; SigInv[2,5] <- 0.00; SigInv[2,6] <- 0.00;
SigInv[3,1] <-  0.00; SigInv[3,2] <-  0.00; SigInv[3,3] <-  1.00; SigInv[3,4] <- -0.40; SigInv[3,5] <- 0.00; SigInv[3,6] <- 0.00;
SigInv[4,1] <-  0.00; SigInv[4,2] <-  0.00; SigInv[4,3] <- -0.40; SigInv[4,4] <-  2.00; SigInv[4,5] <- 0.00; SigInv[4,6] <- 0.00;
SigInv[5,1] <-  0.00; SigInv[5,2] <-  0.00; SigInv[5,3] <-  0.00; SigInv[5,4] <-  0.00; SigInv[5,5] <- 1.00; SigInv[5,6] <- 0.00;
SigInv[6,1] <-  0.00; SigInv[6,2] <-  0.00; SigInv[6,3] <-  0.00; SigInv[6,4] <-  0.00; SigInv[6,5] <- 0.00; SigInv[6,6] <- 1.00;

G <- ( I-A ) 
K.lr  = t(G) %*% SigInv %*% G
dimnames(K.lr) <- list(names,names)

AdjLR  = (K.lr != 0)*1
AdjLR[ row(I)==col(I) ] <- 0
dimnames(AdjLR) <- list(names,names)

# simulate
T <- 1000
y <- matrix(0,T,N,dimnames=list( NULL , names ))
eps <- mvrnorm(T, rep(0,N) , solve(SigInv) )
for( t in 2:T ){ y[t,] = A %*% y[t-1,] + eps[t,] }

# nets procedure
network <- nets(y,p=1,lambda=1)

# plot
plot( network )

#
print( cbind( AdjLR , rep(NA,N) , network$Adj ) )
cat('\n')

print( cbind( round(K.lr,2) , rep(NA,N) , round(network$KLR,2) ) )
cat('\n')
