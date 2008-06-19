\name{lcd.getweights}
\alias{lcd.getweights}

\title{Given a matrix of data, find how many times each occurs}
\description{
  This function takes takes a matrix of binned data and returns a
  matrix containing the distinct observations and a vector of weights}

\usage{ lcd.getweights(x) }

\arguments{
  \item{x}{a data matrix}
  		        }
\details{
    This function is very simple: if we have some data where not all
  the observations are distinct (due to e.g. finite precision
  observation), we can construct a vector of weights for use in
  computing the maximum likelihood estimator.
}
\value{
  \item{x}{The distinct rows of the input matrix}
  \item{w}{Weights: row \eqn{i} appears \eqn{w_i} times}
  }

\author{Madeleine Cule \email{mlc40@cam.ac.uk}
  Robert Gramacy
  Richard Samworth}
\examples{
## simple normal example
set.seed(101)
x <- matrix(rnorm(200),ncol=2)
tmp <- lcd.getweights(x)
out <- lcd.mle(tmp$x,tmp$w)
plot(out,itype="ic")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
