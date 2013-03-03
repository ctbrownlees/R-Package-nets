\name{nets}
\alias{nets}
\title{Network Estimation for Time Series}
\description{`nets` is used to fit sparse long run partial correlation networks for time series using the NETS algorithm as well as other types of sparse networks.}
\usage{ nets( y , type='lrpc' , algorithm='default' , p=1 , lambda=NULL , verbose=FALSE ) }
\arguments{
\item{y}{data, an T x N matrix, each column being a data time series.}
\item{type}{network type: Long Run Partial Correlation (lrpc), Partial Correlation (pc) or Granger (g)}
\item{p}{VAR order (for Granger)}
\item{algorithm}{Optimization algorithm to be used.}
\item{lambda}{shrinkage parameter}
\item{verbose}{extra output messages} }
\details{Details}
\value{
The return value of the nets routine is a network object.
}
\references{
Barigozzi, M. and Brownlees, C. (2012) NETS: Network Estimation for Time Series
}
\author{Christian Brownlees}
\note{
}
\seealso{
}
\examples{
}
\keyword{network}
\keyword{timeseries}
