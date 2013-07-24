\name{getweights}
\alias{getweights}


\title{Find appropriate weights for likelihood calculations}

\description{ This function takes takes a \code{matrix} of (possibly
binned) data and returns a \code{matrix} containing the distinct
observations, and a \code{vector} of weights \eqn{w} as described below.
}

\usage{getweights(x)}

\arguments{
\item{x}{a data \code{matrix}}
}

\details{
  Given an \eqn{n \times d}{n x d} \code{matrix} \eqn{x} of points in
  \eqn{R^d}{R^d}, this function removes duplicated observations, and
  counts the number of times each observation occurs.  This is used to
  compute a \code{vector} \eqn{w} such that \deqn{w_i = \frac{\# \textrm{ of
  times value } i\textrm{ is observed }}{\# \textrm{ of
  observations}}.}{w[i] = (number of times x[i,] is observed)/(total
  number of observations).}

  This function is called by \code{\link{mlelcd}} in order to compute
  the maximum likelihood estimator when the observed data values are not
  distinct.  In this case, the log likelihood function is of the form
  \deqn{\sum_{j=1}^m w_j \log f(X_j),}{w_1 log f(X_1), ..., w_m log f(X_m),}
  where the sum is over distinct observations.

}

\value{ 
  \item{xout}{A \code{matrix} containing the distinct rows of the input
  \code{matrix} \code{x}}
  \item{w}{A real-valued \code{vector} of weights as described above}
}

\author{Madeleine Cule
  
  Robert Gramacy

  Richard Samworth}

\examples{
## simple normal example

x <- matrix(rnorm(200),ncol=2)
tmp <- getweights(x)
lcd <- mlelcd(tmp$x,tmp$w)
plot(lcd,type="ic")
}

\seealso{
  \code{\link{mlelcd}}
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
