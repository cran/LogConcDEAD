\name{rlcd}
\alias{lcd.sample}
\alias{rlcd}

\title{Sample from a log-concave maximum likelihood estimate}

\description{Draws independent samples a log-concave maximum likelihood
  estimate. The estimate should be specified in the form of an object of
  class \code{"LogConcDEAD"}, the result of a call to
  \code{\link{mlelcd}}.
}

\usage{
rlcd(n, lcd)
lcd.sample(lcd, nsample=1) }

\arguments{
\item{n, nsample}{A scalar integer indicating the number of samples required} 
\item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from
  \code{\link{mlelcd}})}
}

\details{
  This function uses a simple rejection sampling scheme to draw
  independent random samples from a log-concave maximum likelihood
  estimator. 

  For examples, see \code{\link{mlelcd}}.
  
\code{lcd.sample} is deprecated, but retained for compatibility with previous versions of this package.
}
\value{
  A numeric \code{matrix} with \code{nsample} rows, each row corresponding to a point
  in \eqn{R^d}{R^d} drawn from the distribution with density defined by \code{lcd}.
}

\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

  Richard Samworth}

\note{Details of the rejection sampling scheme can be found in
  \cite{Cule, Samworth and Stewart (2008)}
  }

  \references{  Cule, M. L., Samworth, R. J. and Stewart, M. I. (2007)
  \emph{Maximum likelihood estimation of a log-concave density}
  Submitted, available from arXiv at \url{http://arxiv.org/abs/0804.3989}
}

\seealso{ \code{\link{mlelcd}}}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{datagen}
\keyword{distribution}
