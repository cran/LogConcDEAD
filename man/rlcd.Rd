\name{rlcd}
\alias{rlcd}

\title{Sample from a log-concave maximum likelihood estimate}

\description{Draws samples from a log-concave maximum likelihood
  estimate. The estimate should be specified in the form of an object of
  class \code{"LogConcDEAD"}, the result of a call to
  \code{\link{mlelcd}}.
}

\usage{
rlcd(n=1, lcd, method=c("Independent","MH"))
}

\arguments{
\item{n}{A scalar integer indicating the number of samples required} 
\item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from
  \code{\link{mlelcd}})}
\item{method}{Indicator of the method used to draw samples, either via independent rejection sampling
  (default choice) or via Metropolis-Hastings} 
}

\details{
  This function by default uses a simple rejection sampling scheme to 
  draw independent random samples from a log-concave maximum likelihood
  estimator. One can also use the Metropolis-Hastings option to draw (dependent)
  samples with a higher acceptance rate.

  For examples, see \code{\link{mlelcd}}.
}

\value{
  A numeric \code{matrix} with \code{nsample} rows, each row corresponding to a point
  in \eqn{R^d}{R^d} drawn from the distribution with density defined by \code{lcd}.
}

\author{
  Yining Chen

  Madeleine Cule

  Robert Gramacy

  Richard Samworth}

\note{Details of the rejection sampling can be found in Appendix B.3 of
  \cite{Cule, Samworth and Stewart (2010)}. Details of the Metropolis-Hastings 
  scheme can be found in \cite{Gopal and Casella (2010)}
  }

\references{
Cule, M. L., Samworth, R. J., and Stewart, M. I. (2010)
 \emph{Maximum likelihood estimation of a multi-dimensional log-concave density} 
 J. Roy. Statist. Soc., Ser. B. (with discussion), 72, 545-600.

Gopal, V. and Casella, G. (2010) 
 \emph{Discussion of Maximum likelihood estimation of a log-concave density by Cule, Samworth and Stewart}  
 J. Roy. Statist. Soc., Ser. B., 72, 580-582.

}

\seealso{ \code{\link{mlelcd}}}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{datagen}
\keyword{distribution}
