\name{interplcd}
\alias{lcd.interp}
\alias{interplcd}

\title{Evaluate the log-concave maximum likelihood estimator of 2-d data on a grid for plotting}
\description{
  Uses the \pkg{akima} package to evaluate the logarithm of the
  log-concave maximum likelihood estimator on a grid for 2-d data, for
  use in \code{\link{plot.LogConcDEAD}}.}

\usage{interplcd(lcd, gridlen=100, span=0.5, ...)
lcd.interp(lcd, gridlen=100, span=0.5, ...) }

\arguments{
  \item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from \code{\link{mlelcd}})}
  \item{gridlen}{A scalar indicating the size of the grid}
  \item{span}{Scalar parameter needed by \code{\link[akima]{interp}}
  (see \code{\link[akima]{interp}} for more details)}
  \item{...}{Other arguments to be passed to \code{\link[akima]{interp}}}
  }
\details{
    Uses the \pkg{\link[akima:akima-package]{akima}} package's
    \code{\link[akima:interp]{interp}} function to evaluate the log-concave maximum likelihood estimator on a grid covering the span of the data.
    
    The output is of a form readily usab;e by
    \code{\link{plot.LogConcDEAD}}, \code{\link[base]{image}},
    \code{\link[base]{contour}}, etc, as illustrated in the examples below.

    For examples, please see \code{\link{mlelcd}}.
    
    \code{lcd.interp} is deprecated, but retained for compatibility with
    previous versions.  }

  \value{
  \item{x}{Vector of \eqn{x}-values of the grid}
  \item{y}{Vector of \eqn{y}-values of the grid}
  \item{z}{A \code{matrix} of the values of the log of the  maximum
  likelihood estimator at points on the grid}
}

\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

Richard Samworth}
\note{
  Needs the \pkg{\link[akima:akima-package]{akima}} package.
}
\seealso{ \code{\link[akima]{interp}} , \code{\link{mlelcd}}}
\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
