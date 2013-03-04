\name{nets}
\alias{nets}
\alias{print.nets}
\alias{plot.nets}
\title{Network Estimation for Time Series}
\description{`nets` is used to fit sparse long run partial correlation networks for time series using the NETS algorithm as well as other types of sparse networks.}
\usage{ nets( y , type='lrpc' , algorithm='default' , p=1 , lambda , verbose=FALSE ) }
\arguments{
\item{y}{data, an T x N matrix, each column being a data time series.}
\item{type}{network type: Long Run Partial Correlation (lrpc), Partial Correlation (pc) or Granger (g)}
\item{p}{VAR order (for Granger)}
\item{algorithm}{Optimization algorithm to be used.}
\item{lambda}{shrinkage parameter(s)}
\item{verbose}{extra output messages} }
\details{The nets procedure can be used to fit long run partial correlation networks, partial correlation and Granger networks.}
\value{The return value of the nets method is a network object.}
\references{
Barigozzi, M. and Brownlees, C. (2012) NETS: Network Estimation for Time Series
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

network <- nets(y,type='g',lambda=1:5,verbose=TRUE)

print( cbind( A , rep(NA,N) , round( network$A[1,,] , 2 ) ) )

}
\keyword{network}
\keyword{multivariate timeseries}
