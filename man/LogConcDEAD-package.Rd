\name{LogConcDEAD-package}
\alias{LogConcDEAD-package}
\alias{LogConcDEAD}

\docType{package}

\title{
Maximum likelihood estimation of a log-concave density in arbitrary dimensions
}

\description{ This package contains a function to compute the maximum
likelihood estimator of a log-concave density in any number of
dimensions using Shor's r-algorithm.

Functions to plot (for 1- and 2-d data) and draw samples from the MLE
are provided.  }

\details{
\code{\link{lcd.mle}} computes the MLE (specified via its value at data points).
Output is a list of class \code{"LogConcDEAD"} which is used for plotting,
function evaluation etc. 

\code{\link{lcd.eval}} evaluates the MLE at a particular point.

\code{\link{lcd.sample}} draws samples from the MLE.

\code{\link{lcd.interp}} interpolates the MLE on a grid, for plotting
purposes.

\code{\link{lcd.marg}} integrates the MLE in 2-d to allow plotting of
estimates of (currently only axis-aligned) marginals.

\code{\link{plot.LogConcDEAD}} produces plots of the MLE, optionally
using the \pkg{rgl} package.


}

\author{
  Madeleine Cule \email{mlc40@cam.ac.uk}
  
  Robert Gramacy \email{bobby@statslab.cam.ac.uk}

  Richard Samworth \email{rjs57@cam.ac.uk}

  Maintainer: Madeleine Cule \email{mlc40@cam.ac.uk}
  }

\note{For one dimensional data, the active set algorithm in
  \pkg{logcondens} is much faster.}

\references{
  Cule, M. L., Samworth, R. J. and Stewart, M. I. (2007)
  \emph{Computing the maximum likelihood estimator of a log-concave
    density}
  In preparation

  Kappel, F. and Kuntsevich, A. V. (2000)
  \emph{An implementation of Shor's r-algorithm}
  Computational Optimization and Applications 15(2) p. 193-205

  \url{http://www.uni-graz.at/imawww/kuntsevich/solvopt/}

  Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T. (1996)
  \emph{The Quickhull algorithm for convex hulls}
  ACM Trans. on Mathematical Software, 22(4) p. 469-483
  \url{http://www.qhull.org}
}

\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ package }

\seealso{
  \code{\link[<logcondens>:<logcondens>-package]{logcondens}}
  \code{\link[<rgl>:<rgl>-package]{rgl}}
}

\examples{

#example with some simple normal data

set.seed(101)
x <- matrix(rnorm(200),ncol=2)
out <- lcd.mle(x)
plot(out,itype="ic")
}
