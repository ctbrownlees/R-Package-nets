\name{nets}
\alias{nets}
\alias{print.nets}
\alias{predict.nets}
\title{Network Estimation for Time Series}
\description{`nets` is used to fit sparse VARs using the NETS algorithm.}
\usage{ nets( y , p=1 , GN=TRUE , CN=TRUE , lambda , alpha.init=NULL , rho.init=NULL , algorithm='activeshooting' , weights='adaptive' , iter.in=100 , iter.out=2 , verbose=FALSE ) }
\arguments{
\item{y}{data, an T x N matrix, each column being a time series.}
\item{p}{VAR order (default 1)}
\item{GN}{Estimate the Autoregressive VAR matrices (default true)}
\item{CN}{Estimate the Concentration matrix of the VAR innovations (default true)}
\item{lambda}{shrinkage parameter either a scalar or a vector of size two}
\item{alpha.init}{initial value of the alpha parameter}
\item{rho.init}{initial value of the rho parameter}
\item{algorithm}{type of algorithm used for the lasso optimization 'shooting' or 'activeshooting' (default 'activeshooting')}
\item{weights}{lasso weights: 'none' or 'adaptive' (default 'adaptive')}
\item{iter.in}{maximum number of in iterations}
\item{iter.out}{maximum number of out iterations}
\item{verbose}{extra output messages} }
\details{The nets procedure estimates sparse VAR models by LASSO.}
\value{
	The nets procedure returns an object of class "nets"

	An object of class "nets" is a list containing at least the following components:

	A.hat: an N x N x P array containing the estimated autoregressive matrices
  
  C.hat: an N x N matrix containing the estimated concentration matrix

  The function "print" is used to obtain and print a summary of the estimation results.
}
\references{
Barigozzi, M. and Brownlees, C. (2016) NETS: Network Estimation for Time Series, Working Paper
}
\author{Christian Brownlees}
\examples{
1+1
}
\keyword{network}
\keyword{multivariate timeseries}
