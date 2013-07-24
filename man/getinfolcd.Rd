\name{getinfolcd}
\Rdversion{1.1}
\alias{getinfolcd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct an object of class LogConcDEAD
}
\description{
  A function to construct an object of class \code{LogConcDEAD} from a
  dataset (given as a matrix) and the value of the log maximum
  likelihood estimator at datapoints.
}
\usage{
getinfolcd(x, y, w = rep(1/length(y), length(y)), chtol = 10^-6, 
  MinSigma = NA, NumberOfEvaluations = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Data in \eqn{R^d}, in the form of an \eqn{ n \times d}{n x d}
    numeric \code{matrix}
}
\item{y}{
  Value of log of maximum likelihood estimator at data points
}
\item{w}{
  Vector of weights \eqn{w_i} such that the computed estimator
  maximizes \deqn{\sum_{i=1}^n w_i \log f(x_i)}{w[1] log f(x[1,]) + ... +
  w[n] log f([x,n])} subject to the restriction that \eqn{f} is
  log-concave. The default is  \eqn{\frac{1}{n}}{1/n} for all \eqn{i},
  which corresponds to i.i.d. observations.   
}
\item{chtol}{
  Tolerance for computation of convex hull. Altering this is not recommended.
}
\item{MinSigma}{
  Real-valued scalar giving minimum value of the objective function
}
  \item{NumberOfEvaluations}{Vector containing the number of steps, number of function
    evaluations, and number of subgradient evaluations.  If the \pkg{SolvOpt}
    algorithm fails, the first component will be an error code \eqn{(<0)}
  }
 } 
\details{
  This function is used in \code{\link{mlelcd}}
}
\value{ An object of class \code{"LogConcDEAD"}, with the following
  components:
  \item{x}{Data copied from input (may be reordered)}
  \item{w}{weights copied from input (may be reordered)}
  \item{logMLE}{\code{vector} of
    the log of the maximum likelihood estimate, evaluated at the observation points}
  \item{NumberOfEvaluations}{Vector containing the number of steps, number of function
    evaluations, and number of subgradient evaluations.  If the \pkg{SolvOpt}
    algorithm fails, the first component will be an error code \eqn{(<0)}.}
  \item{MinSigma}{Real-valued scalar giving minimum value of the objective function}
  \item{b}{\code{matrix} (see Details)}
  \item{beta}{\code{vector} (see Details)}
  \item{triang}{\code{matrix} containing final triangulation of the convex hull of the data}
  \item{verts}{\code{matrix} containing details of triangulation for use in \code{\link{dlcd}}}
  \item{vertsoffset}{\code{matrix} containing details of triangulation for use in \code{\link{dlcd}}}
\item{chull}{Vector containing vertices of faces of the convex hull of
  the data}
\item{outnorm}{\code{matrix} where each row is an outward
  pointing normal vectors for the faces of the convex hull of the
  data. The number of vectors depends on the number of faces of the
  convex hull.}
\item{outoffset}{\code{matrix} where each row is a point on a face of
  the convex hull of the data. The number of vectors depends on the
  number of faces of the convex hull.}
}

\author{Madeleine Cule
  
  Robert B. Gramacy

  Richard Samworth

  Yining Chen
}


\seealso{
  \code{\link{mlelcd}}
}

\keyword{ nonparametric }
\keyword{ multivariate }
