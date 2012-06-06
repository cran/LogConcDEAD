\name{interplcd}
\alias{interplcd}


\title{Evaluate the log-concave maximum likelihood estimator of 2-d data
  on a grid for plotting}

\description{ Evaluates the logarithm of the
log-concave maximum likelihood estimator on a grid for 2-d data, for use
in \code{\link{plot.LogConcDEAD}}.}

\usage{interplcd(lcd, gridlen=100 )}

\arguments{
  \item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from \code{\link{mlelcd}})}
  \item{gridlen}{A scalar indicating the size of the grid}

}
\details{
  Interpolates the MLE over a grid.
  
    The output is of a form readily usable by
    \code{\link{plot.LogConcDEAD}}, \code{\link[graphics]{image}},
    \code{\link[graphics]{contour}}, etc, as illustrated in the examples below.

    For examples, please see \code{\link{mlelcd}}.
    
}

  \value{
  \item{x}{Vector of \eqn{x}-values of the grid}
  \item{y}{Vector of \eqn{y}-values of the grid}
  \item{z}{A \code{matrix} of the values of the log of the  maximum
  likelihood estimator at points on the grid}
}

\author{Madeleine Cule

  Robert Gramacy

Richard Samworth}
\seealso{ \code{\link{mlelcd}}}
\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
