\name{lcd.interp}
\alias{lcd.interp}

\title{Evaluate the log-concave MLE of 2-d data on a grid for plotting}
\description{
  A function (using the \pkg{akima} package) to evaluate the logarithm of the
  log-concave MLE on a grid for 2-d data, for use in \code{\link{plot.LogConcDEAD}}.}

\usage{ lcd.interp(out, gridlen=100, span=0.5, ...) }

\arguments{
  \item{out}{Output from lcd.mle (of class LogConcDEAD)}
  \item{gridlen}{Size of grid}
  \item{span}{Parameter needed by \code{interp}}
    \item{...}{Other arguments to be passed to \code{interp}}
    }
\details{
    Uses the \pkg{akima} package's \code{interp} function to linearly interpolate.
}
\value{
  \item{x}{x-values of grid}
  \item{y}{y-values of grid}
  \item{z}{Log MLE on a grid}
 
}

\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

Richard Samworth}
\note{
  Needs the \pkg{akima} package.
}
\seealso{ \code{\link[<akima>:<akima>-package]{<interp>}}
          \code{\link[<rgl>:<rgl>-package]{<rgl>}}}
\examples{
#simple normal data, a couple of plots

x <- matrix(rnorm(200),ncol=2)
out <- lcd.mle(x)
g <- lcd.interp(out)
par(mfrow=c(2,2))
plot(out, g=g, itype="c", main="Density")
plot(out, g=g, itype="c", uselog=TRUE, main="LogDensity")
plot(out, g=g, itype="i", main="Density")
plot(out, g=g, itype="i", uselog=TRUE, main="LogDensity")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
