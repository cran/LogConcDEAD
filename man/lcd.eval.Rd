\name{lcd.eval}
\alias{lcd.eval}

\title{Evaluate the log-concave MLE at a point}
\description{
 Function to evaluate the log-concave MLE at a point}
\usage{
lcd.eval(out, po, uselog=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{out}{Object of class \code{LogConcDEAD} (output from \code{\link{lcd.mle}})}
  \item{po}{Point (or matrix of points) at which MLE should be evaluated}
  \item{uselog}{Logical for whether the estimator should be calculated
    on the log scale}
}
\details{
  The function is evaluated using the special structure of the MLE.  The
  appropriate triangulation is returned by \code{\link{lcd.mle}} along with the
  equation of the hyperplanes defining the surface of the log MLE.
}
\value{
 The MLE (or log MLE) evaluated at the point.}
\author{Madeleine Cule \email{mlc40@cam.ac.uk}
  
  Robert Gramacy
  
  Richard Samworth}

\seealso{ \code{\link{lcd.mle}}}
\examples{
set.seed(101)
x <- matrix(rnorm(200),ncol=2)
out <- lcd.mle(x)
point <- c(0,0)
lcd.eval(out,point,uselog=FALSE)
}

\keyword{multivariate}
\keyword{nonparametric}
