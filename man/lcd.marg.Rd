\name{lcd.marg}
\alias{lcd.marg}
\alias{lcd.marg2}

\title{Find axis-aligned marginals of the MLE for 2-d data}
\description{
  This function numerically integrates the log MLE of 2-d data to
  produce axis-aligned marginals for use in \code{\link{plot.LogConcDEAD}}.}

\usage{ lcd.marg (out, gridlen=100)

lcd.marg2 (out, gridlen=100) }

\arguments{
  \item{out}{Output from lcd.mle (of class LogConcDEAD)}
  \item{gridlen}{Size of grid}
  }
\details{
    Does a (currently rather crude) approximation to the integral of the
    log MLE, for the purpose of plotting axis-aligned marginals. The
    \code{integrate} routine is used to perform numerical integration.

    Currently use \code{lcd.marg} for \eqn{X_1} and \code{lcd.marg2} for
    \eqn{X_2}{X_2}.
}

\value{
  \item{x}{x-values of grid (for \eqn{X_1} axis)}
  \item{y}{y-values of grid (for \eqn{X_2} axis)}
  \item{marg}{Integrated MLE}
}


\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

  Richard Samworth}


\examples{
#simple normal data, a couple of plots

x <- matrix(rnorm(200),ncol=2)
out <- lcd.mle(x)
m1 <- lcd.marg(out)
m2 <- lcd.marg2(out)
par(mfrow=c(1,2))
plot(out,marg=1,g.marg=m1)
plot(out,marg=2,g.marg=m2)
}

\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
