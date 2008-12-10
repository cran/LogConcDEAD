\name{dmarglcd}
\alias{lcd.marg.eval}
\alias{dmarglcd}
\title{Evaluate the marginal of multivariate log-concave maximum likelihood estimators at a point}
\description{
  Integrates the log-concave maximum likelihood estimator of multivariate data to
  evaluate the marginal density at a point.}

\usage{dmarglcd(x=0, lcd, marg=1)
lcd.marg.eval (lcd, point=0, marg=1)
 }
\arguments{
  \item{x, point}{Point (or \code{vector} of points) at which the marginal density is to be evaluated}
  \item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from \code{\link{mlelcd}})}
  \item{marg}{Which margin is required?}
    }

\details{
      
Given a multivariate log-concave maximum likelihood estimator in the
form of an object of class \code{"LogConcDEAD"}, a margin \code{marg},
and a real-valued point \code{x}, this function evaluates the estimated
marginal density \eqn{\hat{f}_{n,\tiny\texttt{marg}}(x)}{f_n,marg (x)},
as obtained by integrating over all the other dimensions.

For examples, see \code{\link{mlelcd}}.

\code{lcd.marg.eval} is deprecated, but retained for compatibility with previous versions.
   }

\value{
  A \code{vector} containing the values of the marginal density \eqn{\hat{f}_{n,
  \tiny\texttt{marg}}}{f_n,marg} at the points \code{x}.
}

\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

  Richard Samworth}

\seealso{\code{\link{mlelcd}}}

\keyword{dplot}
\keyword{multivariate}
\keyword{nonparametric}
