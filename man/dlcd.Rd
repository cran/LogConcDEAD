\name{dlcd}
\alias{dlcd}
\alias{lcd.eval}

\title{Evaluation of a log-concave maximum likelihood estimator at a
  point}

\description{ This function evaluates the density function of a
log-concave maximum likelihood estimator at a point or points.}

\usage{dlcd(x,lcd, uselog=FALSE) 
lcd.eval(lcd, po, uselog=FALSE) }

\arguments{ 
\item{x, po}{ Point (or \code{matrix} of points) at which the maximum
  likelihood estimator should be evaluated}
\item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from
  \code{\link{mlelcd}})} 
\item{uselog}{Scalar \code{logical}: should the estimator should be calculated on the
  log scale?} 
}

\details{
  A log-concave maximum likelihood estimate
  \eqn{\hat{f}_n}{f_n} is satisfies \eqn{\log \hat{f}_n = \bar{h}_y}{log
    f_n = h_y} for some \eqn{y \in R^n}{y in R^n}, where
  \deqn{\bar{h}_y(x) = \inf \lbrace h(x) \colon h \textrm{ concave }, h(x_i)
  \geq y_i \textrm{ for } i = 1, \ldots, n \rbrace.}{h_y(x)
  = inf{h(x): h concave, h(x_i) >= y_i for i = 1, ..., n}.}

  Functions of this form may equivalently be specified by dividing
  \eqn{C_n}, the convex hull of the data into simplices \eqn{C_j} for
  \eqn{j \in J} (triangles in 2d, tetrahedra in 3d etc), and setting
  \deqn{f(x) = \exp\{b_j^T x - \beta_j\}}{f(x) = exp{b_j^T x - beta_j}}
  for \eqn{x \in C_j}{x in C_j}, and \eqn{f(x) = 0} for \eqn{x \notin
  C_n}{x not in C_n}.  The estimated density is zero outside the convex
  hull of the data.
  
  The estimate may therefore be evaluated by finding the appropriate
  simplex \eqn{C_j}, then evaluating \eqn{\exp\{b_j^T x -
  \beta_j\}}{f(x) = exp{b_j^T x - beta_j}} (if \eqn{x \notin C_n}{x not
  in C_n}, set \eqn{f(x) = 0}).

For examples, see \code{\link{mlelcd}}.

\code{lcd.eval} is deprecated, but retained for compatibility with
previous versions.  }

\value{ A \code{vector} of  maximum likelihood estimate (or log
  maximum likelihood estimate) values, as evaluated at the points \code{x}.}

\author{
  Madeleine Cule \email{mlc40@cam.ac.uk}
  
  Robert Gramacy
  
  Richard Samworth
}

\seealso{ \code{\link{mlelcd}}}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{distribution}
