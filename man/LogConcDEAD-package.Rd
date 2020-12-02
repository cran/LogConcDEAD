\name{LogConcDEAD-package}
\alias{LogConcDEAD-package}
\alias{LogConcDEAD}

\docType{package}

\title{
  Computes a log-concave (maximum likelihood) estimator for i.i.d. data in any number of dimensions
}

\description{ This package contains a function to compute the maximum
likelihood estimator of a log-concave density in any number of
dimensions using Shor's \eqn{r}-algorithm.

Functions to plot (for 1- and 2-d data), evaluate and draw samples from the maximum likelihood estimator
are provided.  }

\details{
  This package contains a selection of functions for maximum likelihood
  estimation under the constraint of log-concavity.
  
\code{\link{mlelcd}} computes the maximum likelihood estimator
(specified via its value at data points).  Output is a list of class
\code{"LogConcDEAD"} which is used as input to various auxiliary functions.

\code{\link{hatA}} calculates the difference between the sample covariance and the fitted covariance.

\code{\link{dlcd}} evaluates the estimated density at a particular point.

\code{\link{dslcd}} evaluates the smoothed version of estimated density at a particular point.

\code{\link{rlcd}} draws samples from the estimated density.

\code{\link{rslcd}} draws samples from the smoothed version of estimated density.

\code{\link{interplcd}} interpolates the estimated density on
a grid for plotting purposes.

\code{\link{dmarglcd}} evaluates the estimated marginal density by
integrating the estimated density over an appropriate subspace.

\code{\link{interpmarglcd}} evaluates a marginal density estimate at
equally spaced points along the axis for plotting purposes.  This is
done by integrating the estimated density over an appropriate subspace.

\code{\link{plot.LogConcDEAD}} produces plots of the maximum likelihood
estimator, optionally using the \pkg{rgl} package.

\code{\link[=print.LogConcDEAD]{print}} and
\code{\link[=summary.LogConcDEAD]{summary}} methods are also available.
}

\author{
  Yining Chen (maintainer) \email{y.chen101@lse.ac.uk}

  Madeleine Cule 
  
  Robert Gramacy 

  Richard Samworth
}

\note{The authors gratefully acknowledge the assistance of Lutz Duembgen
  at the University of Bern for his insight into the objective function
  in \link{mlelcd}.

  For one dimensional data, the active set algorithm in
  \code{\link[logcondens:logcondens-package]{logcondens}} is much faster.}

\references{
Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T. (1996)
 \emph{The Quickhull algorithm for convex hulls} 
 ACM Trans. on Mathematical Software, 22(4) p.469-483
 \url{http://www.qhull.org}

Chen, Y. and Samworth, R. J. (2013) 
 \emph{Smoothed log-concave maximum likelihood estimation with applications} 
 Statist. Sinica, 23, 1373-1398. \url{https://arxiv.org/abs/1102.1191v4}

Cule, M. L. and D\"umbgen, L. (2008)
 \emph{On an auxiliary function for log-density estimation}, 
 University of Bern technical report. \url{https://arxiv.org/abs/0807.4719}

Cule, M. L., Samworth, R. J., and Stewart, M. I. (2010)
 \emph{Maximum likelihood estimation of a multi-dimensional log-concave density} 
 J. Roy. Statist. Soc., Ser. B. (with discussion), 72, 545-600.

Gopal, V. and Casella, G. (2010) 
 \emph{Discussion of Maximum likelihood estimation of a log-concave density by Cule, Samworth and Stewart}  
 J. Roy. Statist. Soc., Ser. B., 72, 580-582.

Grundmann, A. and Moeller, M. (1978)
 \emph{Invariant Integration Formulas for the N-Simplex by Combinatorial Methods} 
 SIAM Journal on Numerical Analysis, Volume 15, Number 2, 282-290.

Kappel, F. and Kuntsevich, A. V. (2000)
 \emph{An implementation of Shor's r-algorithm} Computational 
 Optimization and Applications, Volume 15, Issue 2, 193-205.

Shor, N. Z. (1985)
 \emph{Minimization methods for nondifferentiable functions}
 Springer-Verlag
}

\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ package }
\keyword{ smoothing }

\seealso{
  \code{\link[logcondens:logcondens-package]{logcondens}},
  \code{\link[rgl:rgl-package]{rgl}}
}

\examples{
## Some simple normal data, and a few plots

x <- matrix(rnorm(200),ncol=2)
lcd <- mlelcd(x)
g <- interplcd(lcd)
par(mfrow=c(2,2), ask=TRUE)
plot(lcd, g=g, type="c")
plot(lcd, g=g, type="c", uselog=TRUE)
plot(lcd, g=g, type="i")
plot(lcd, g=g, type="i", uselog=TRUE)

## Some plots of marginal estimates
par(mfrow=c(1,1))
g.marg1 <- interpmarglcd(lcd, marg=1)
g.marg2 <- interpmarglcd(lcd, marg=2)
plot(lcd, marg=1, g.marg=g.marg1)
plot(lcd, marg=2, g.marg=g.marg2) 

## generate some points from the fitted density
generated <- rlcd(100, lcd)
genmean <- colMeans(generated)

## evaluate the fitted density
mypoint <- c(0, 0)
dlcd(mypoint, lcd, uselog=FALSE)
mypoint <- c(10, 0)
dlcd(mypoint, lcd, uselog=FALSE)

## evaluate the marginal density
dmarglcd(0, lcd, marg=1)
dmarglcd(1, lcd, marg=2)
}
