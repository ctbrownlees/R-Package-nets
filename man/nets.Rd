\name{nets}
\alias{nets}
\alias{print.nets}
\title{Network Estimation for Time Series}
\description{`nets` is used to fit sparse VARs using the NETS algorithm.}
\usage{ nets( y , p=1 , GN=TRUE , CN=TRUE , lambda , alpha.init=NULL , rho.init=NULL , algorithm='activeshooting' , weights='adaptive' , maxiter=100 , verbose=FALSE ) }
\arguments{
\item{y}{data, an T x N matrix, each column being a time series.}
\item{p}{VAR order (default 1)}
\item{GN}{Estimate the Autoregressive VAR matrices (default true)}
\item{CN}{Estimate the Concentration matrix of the VAR innovations (default true)}
\item{lambda}{shrinkage parameter}
\item{alpha.init}{initial value of the alpha parameter}
\item{rho.init}{initial value of the rho parameter}
\item{weights}{lasso weights: 'none' or 'adaptive' (default 'adaptive')}
\item{maxiter}{maximum number of iterations}
\item{verbose}{extra output messages} }
\details{The nets procedure estimates sparse VARs. The routine can be used to estimate long run partial correlation networks, partial correlation networks and Granger networks.}
\value{
	The ‘nets’ procedure returns an object of class ‘"nets"’

	The functions ‘summary’, ‘plot’ and ‘print’ are used to obtain and print a
	summary of the results.

	An object of class ‘"nets"’ is a list containing at least the following components:

	A.hat: an N x N x P array of VAR matrix parameters
  
  C.hat:
}
\references{
Barigozzi, M. and Brownlees, C. (2015) NETS: Network Estimation for Time Series
}
\author{Christian Brownlees}
\examples{

N   <- 3
T   <- 500
 
# A
A      <- matrix( 0 , N , N )
A[1,1] <- 0.71; A[1,2] <- 0.00; A[1,3] <- 0.00; 
A[2,1] <- 0.00; A[2,2] <- 0.63; A[2,3] <- 0.00; 
A[3,1] <- 0.00; A[3,2] <- 0.00; A[3,3] <- 0.10; 

# Simulate Process
y <- matrix(0,T,N)
eps <- matrix( rnorm(T*N,0,1) , T , N )
for( t in 2:T ){ 
	for( i in 1:N ) {
		y[t,i] = sum( A[i,] * y[t-1,] ) + eps[t,i] 
	}
}

network <- nets(y,p=1,lambda=1,verbose=TRUE)

#network

#plot( network )

}
\keyword{network}
\keyword{multivariate timeseries}
