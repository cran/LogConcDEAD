\name{print.LogConcDEAD}

\alias{print.LogConcDEAD}
\alias{summary.LogConcDEAD}

\title{Summarizing log-concave maximum likelihood estimator}

\description{Generic \code{print} and \code{summary} method for objects
  of class \code{"LogConcDEAD"}}
  
  \usage{
\method{print}{LogConcDEAD}(x, ...)
\method{summary}{LogConcDEAD}(object, ...)
    }

\arguments{
  \item{x}{Object of class \code{"LogConcDEAD"} (typically output from
    \code{\link{mlelcd}}), as required by \code{print} }
  \item{object}{Object of class \code{"LogConcDEAD"} (typically output
  from \code{\link{mlelcd}}), as required by \code{summary} }
  \item{...}{Other arguments passed to \code{print} or \code{summary}}
}

\details{
  \code{print} and \code{summary} currently perform the same function.

  If there has been an error computing the maximum likelihood
  estimator, an error message is printed.

  Otherwise, the value of the log maximum likelihood estimator at
  observation points is printed.  The number of interations required by
  the subgradient and the number of function evaluations are also
  printed.}

\author{
  Madeleine Cule
  
  Robert B. Gramacy
  
  Richard Samworth
}

\seealso{\code{\link{mlelcd}}}

\keyword{nonparametric}
\keyword{multivariate}
