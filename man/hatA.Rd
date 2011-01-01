\name{hatA}
\alias{hatA}

\title{Compute the smoothing matrix of the smoothed log-concave maximum likelihood estimator}

\description{ This function computes the matrix \eqn{\hat{A}} of the smoothed log-concave maximum 
  likelihood estimator}

\usage{hatA(lcd)}

\arguments{ 
\item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from
  \code{\link{mlelcd}})} 
}

\details{
  This function evaluates the the matrix \eqn{\hat{A}} of the smoothed log-concave maximum 
  likelihood estimator, which is positive definite, and equals the difference between the 
  sample covariance matrix and the covariance matrix of the corresponding log-concave maximum
  likelihood density estimator.
}

\value{ A \code{matrix} equals \eqn{\hat{A}} of the smoothed log-concave maximum 
  likelihood estimator}

\author{
  Yining Chen

  Madeleine Cule
  
  Robert Gramacy
  
  Richard Samworth
}

\seealso{ \code{\link{cov.LogConcDEAD}}}
\keyword{ nonparametric }
\keyword{ multivariate }
