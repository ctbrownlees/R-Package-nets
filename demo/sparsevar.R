
require(nets)

# Problem Parameters
N   <- 6
T   <- 500
 
# A
A      <- matrix( 0 , N , N )
A[1,1] <- 0.71; A[1,2] <- 0.00; A[1,3] <- 0.00; A[1,4] <- 0.00; A[1,5] <- 0.00; A[1,6] <- 0.15;
A[2,1] <- 0.00; A[2,2] <- 0.63; A[2,3] <- 0.00; A[2,4] <- 0.30; A[2,5] <- 0.00; A[2,6] <- 0.00;
A[3,1] <- 0.00; A[3,2] <- 0.00; A[3,3] <- 0.10; A[3,4] <- 0.00; A[3,5] <- 0.00; A[3,6] <- 0.00;
A[4,1] <- 0.00; A[4,2] <- 0.00; A[4,3] <- 0.10; A[4,4] <- 0.22; A[4,5] <- 0.25; A[4,6] <- 0.00;
A[5,1] <- 0.10; A[5,2] <- 0.00; A[5,3] <- 0.00; A[5,4] <- 0.00; A[5,5] <- 0.34; A[5,6] <- 0.00;
A[6,1] <- 0.00; A[6,2] <- 0.00; A[6,3] <- 0.00; A[6,4] <- 0.00; A[6,5] <- 0.00; A[6,6] <- 0.42;

# Simulate Process
y <- matrix(0,T,N)
eps <- matrix( rnorm( T*N , 0 , 1 ) , T , N )

for( t in 2:T ){ y[t,] = A %*% y[t-1,] + eps[t,] }

results <- nets.var.sfit(y,lambda.range=c(seq(0.05,20,0.25)),verbose=TRUE)

print( cbind( A , rep(NA,N) , round( results$A[1,,] , 2 ) ) )

mse <- sqrt( sum( (A - results$A[1,,])^2 ) )
tp1 <- mean( (results$A[1,,])[ A==0 ] != 0 )
pow <- mean( (results$A[1,,])[ A!=0 ] != 0 )

cat( 'MSE' , mse , 'TYPE1' , tp1 , 'POW' , pow , 'lambda' , results$lambda , '\n')

