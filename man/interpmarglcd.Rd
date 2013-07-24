\name{interpmarglcd}
\alias{interpmarglcd}


\title{Finds marginals of multivariate logconcave maximum likelihood estimators by integrating}
\description{
  Integrates the maximum likelihood estimator of multivariate data over an appropriate subspace to
  produce axis-aligned marginals for use in \code{\link{plot.LogConcDEAD}}.}

\usage{interpmarglcd(lcd, marg=1, gridlen=100)}

\arguments{
  \item{lcd}{Output from \code{\link{mlelcd}} (of class \code{"LogConcDEAD")}}
  \item{marg}{An (integer) scalar indicating which margin is required}
  \item{gridlen}{An (integer) scalar indicating the size of the grid}
  }

  \details{ Given a multivariate log-concave maximum likelihood
    estimator in the form of an object of class \code{"LogConcDEAD"} and
    a margin \code{marg}, this function will compute the marginal
    density estimate \eqn{\hat{f}_{n, \textrm{\tiny{marg}}}}{f_n,marg}.
    The estimate is evaluated at \code{gridlen} equally spaced points in
    the range where the density estimate is nonzero. These points are
    given in the \code{vector} \code{xo}.

    \eqn{\hat{f}_{n, \textrm{\tiny{marg}}}}{f_n,marg} is evaluated by
    integrating the log-concave maximum likelihood estimator
    \eqn{\hat{f}_n}{f_n} over the other components.  The marginal
    density is zero outside the range of \code{xo}.

    For examples, see \code{\link{mlelcd}}.
    
}

\value{
  \item{xo}{Vector of values at which the marginal density is estimate
  is computed. }
  \item{marg}{Vector of values of the integrated maximum likelihood
  estimator at the locations \code{xo}}
}

\seealso{\code{\link{dmarglcd}}, \code{\link{mlelcd}}}
\author{Madeleine Cule

  Robert Gramacy

  Richard Samworth}

\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
