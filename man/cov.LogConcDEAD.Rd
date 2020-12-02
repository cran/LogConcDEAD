\name{cov.LogConcDEAD}
\alias{cov.LogConcDEAD}

\title{Compute the covariance matrix of a log-concave maximum likelihood estimator}

\description{ This function computes the covariance matrix of a
log-concave maximum likelihood estimator.}

\usage{cov.LogConcDEAD(lcd)}

\arguments{ 
\item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from
  \code{\link{mlelcd}})} 
}

\details{
  This function evaluates the covariance matrix of a given log-concave maximum likelihood
  estimator using the second order partial derivatives of the auxiliary function studied 
  in \cite{Cule, M. L. and D\"umbgen, L. (2008)}.

  For examples, see \code{\link{mlelcd}}.
}

\value{ A \code{matrix} equals the covariance matrix of the 
  log-concave maximum likelihood density estimator.}

\references{
Cule, M. L. and D\"umbgen, L. (2008) On an auxiliary function for
log-density estimation, University of Bern technical report. \url{https://arxiv.org/abs/0807.4719}
}

\author{
  Yining Chen

  Madeleine Cule
  
  Robert Gramacy
  
  Richard Samworth
}

\seealso{ \code{\link{hatA}}}
\keyword{ nonparametric }
\keyword{ multivariate }
