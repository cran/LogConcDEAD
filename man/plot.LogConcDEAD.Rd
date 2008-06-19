\name{plot.LogConcDEAD}
\alias{plot.LogConcDEAD}

\title{Plot function for log-concave MLE}
\description{
 This plot function currently only plots 1- and 2-d data.  }
\usage{
\method{plot}{LogConcDEAD}(x, uselog=FALSE, method="akima", itype="ic", addp=TRUE, gridlen=100,
g=NULL, marg=NULL, g.marg=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{LogConcDEAD} (e.g. as output from lcd.mle)}
  \item{uselog}{logical indicating whether the log or true scale should
    be plotted}
  \item{method}{The \pkg{akima} package is used for linear interpolation}
  \item{itype}{Plot type: \code{p} perspectice, \code{c} contour, \code{i} image, \code{ic}
    image and contour, \code{r} using \pkg{rgl} (the best!)}
  \item{addp}{logical: should the data points be plotted on the surface?}
  \item{gridlen}{Number of points at which the MLE is evaluated in each direction}
  \item{g}{If many plots are required, it may be quicker to compute the
    grid using \code{\link{lcd.interp}(out)} and pass the result to \code{plot}}
  \item{marg}{axis-aligned marginals can be plotted (set to 1 or 2)}
  \item{g.marg}{may be quicker to compute marginal using \code{\link{lcd.marg}}
     and pass to \code{plot}}
\item{...}{Other \code{plot} arguments}
}
\details{
  This function makes use of the \pkg{akima} package.

  For one-dimensional data, the log MLE is piecewise linear with knots
  at observation points.

  The same is (roughly speaking) true in higher dimensions.  The MLE is
evaluated on a grid
  using the \code{interp} function from \pkg{akima}.

  To come: plotting of marginals for higher-dimensional data.
}


\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert B. Gramacy

  Richard Samworth}

\seealso{\code{\link[akima]{interp}}}

\examples{

## some simple normal data

set.seed(101)
x <- matrix(rnorm(200), ncol=2)
out <- lcd.mle(x)
g <- lcd.interp(out)
par(mfrow=c(2,2))
plot(out, g=g, itype="c")
plot(out, g=g, itype="c", uselog=TRUE)
plot(out, g=g, itype="i")
plot(out, g=g, itype="i", uselog=TRUE)
}


\keyword{ dplot }
\keyword{ hplot }
\keyword{ nonparametric }
\keyword{ multivariate }
