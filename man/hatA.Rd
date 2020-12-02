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
  sample covariance matrix and the covariance matrix of the fitted log-concave maximum
  likelihood density estimator.

  For examples, see \code{\link{mlelcd}}
}

\value{ A \code{matrix} equals \eqn{\hat{A}} of the smoothed log-concave maximum 
  likelihood estimator}

\author{
  Yining Chen

  Madeleine Cule
  
  Robert Gramacy
  
  Richard Samworth
}

\note{Details of the computational aspects can be found in \cite{Chen and Samworth (2011)}.}

\references{
 Chen, Y. and Samworth, R. J. (2013)
 \emph{Smoothed log-concave maximum likelihood estimation with applications} 
 Statist. Sinica, 23, 1373-1398. \url{https://arxiv.org/abs/1102.1191v4}
}


\seealso{ \code{\link{cov.LogConcDEAD}}}
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ smoothing }

