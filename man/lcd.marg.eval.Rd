\name{lcd.marg.eval}
\alias{lcd.marg.eval}

\title{Evaluates the marginal of multivariate logconcave MLEs at a point by integrating}
\description{
  This function integrates the MLE of multivariate data to
  evaluate the marginal density at a point.}

\usage{ lcd.marg.eval (out, point=0, marg=1)
}
\arguments{
  \item{out}{Output from lcd.mle (of class \code{LogConcDEAD)}}
  \item{point}{Point at which the density is to be evaluated}
  \item{marg}{Which margin is required?}
    }

\details{
Given a multivariate log-concave MLE in the form of an object of class \code{LogConcDEAD} and a margin \code{marg}, and a real-valued point \code{point}, this function will integrate the MLE over the relevant subspaces.  

The \pkg{R} code for computing the function \code{J} was provided by Lutz D\"umbgen.	
   }

\value{
  The value of the marginal at that point.
}


\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

  Richard Samworth}


\examples{
##simple normal data

set.seed(101)
x <- matrix(rnorm(200),ncol=2)
out <- lcd.mle(x)
lcd.marg.eval(out, 0, 1)
lcd.marg.eval(out, 1, 2)
}

\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
