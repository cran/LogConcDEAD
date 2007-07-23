\name{lcd.sample}
\alias{lcd.sample}

\title{Sample from the MLE of a log-concave density using rejection sampling}
\description{
 Draws independent samples from the MLE of a log-concave density.
}
\usage{
lcd.sample(out, nsample=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{Object of class \code{LogConcDEAD} (output from \code{\link{lcd.mle}})}
  \item{nsample}{Number of samples required}
}
\details{
  A simple rejection sampling scheme is used to draw from the MLE.  For
  discussion of the structure of the MLE, see \code{\link{lcd.mle}}.

  1. A simplex is selected
 
  2. Rejection sampling is used to sample from the MLE restricted to
  that simplex
}
\value{
  A numeric matrix with nsample rows, each row corresponding to a point
  in \eqn{R^d}{R^d}.
}

\author{Madeleine Cule \email{mlc40@cam.ac.uk}

  Robert Gramacy

  Richard Samworth}

\seealso{ \code{\link{lcd.mle}}}
\examples{
#simple normal data

set.seed(101)
x <- matrix(rnorm(200),ncol=2)
out <- lcd.mle(x)
lcd.sample(out,nsample=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{nonparametric}
