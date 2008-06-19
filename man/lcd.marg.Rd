\name{lcd.marg}
\alias{lcd.marg}

\title{Finds marginals of multivariate logconcave MLEs by integrating}
\description{
  This function integrates the MLE of multivariate data to
  produce axis-aligned marginals for use in \code{\link{plot.LogConcDEAD}}.}

\usage{ lcd.marg (out, marg=1, gridlen=100)
}
\arguments{
  \item{out}{Output from lcd.mle (of class \code{LogConcDEAD)}}
  \item{marg}{Which margin is required?}
  \item{gridlen}{Size of grid}
  }

\details{
Given a multivariate log-concave MLE in the form of an object of class \code{LogConcDEAD} and a margin \code{marg}, this will integrate the MLE over the relevant subspaces along a line with \code{gridlen} points.  

The \pkg{R} code for computing the function \code{J} was provided by Lutz D\"umbgen.	
   }

\value{
  \item{xo}{x-values}
  \item{marg}{Integrated MLE}
}


\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

  Richard Samworth}


\examples{
#simple normal data, a couple of plots

set.seed(101)
x <- matrix(rnorm(200),ncol=2)
out <- lcd.mle(x)
m1 <- lcd.marg(out,marg=1)
m2 <- lcd.marg(out,marg=2)
par(mfrow=c(1,2))
plot(out,marg=1,g.marg=m1)
plot(out,marg=2,g.marg=m2)
}

\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
