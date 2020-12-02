\name{dslcd}
\alias{dslcd}

\title{Evaluation of a smoothed log-concave maximum likelihood estimator at given
  points}

\description{ This function evaluates the density function of a
smoothed log-concave maximum likelihood estimator at a point or points.}

\usage{dslcd(x, lcd, A=hatA(lcd))}

\arguments{ 
\item{x}{ Point (or \code{matrix} of points) at which the smoothed log-concave 
  maximum likelihood estimator should be evaluated}
\item{lcd}{Object of class \code{"LogConcDEAD"} (typically output from
  \code{\link{mlelcd}})} 
\item{A}{A positive definite \code{matrix} that determines the degree of smoothing, 
  typically taken as the output of \code{hatA(lcd)}}
}

\details{
  The smoothed log-concave maximum likelihood estimator is a fully automatic
  nonparametric density estimator, obtained as a canonical smoothing of the 
  log-concave maximum likelihood estimator. More precisely, it equals the convolution
  \eqn{ \hat{f} * \phi_{d,\hat{A}}}, where \eqn{\phi_{d,\hat{A}}} is the density function of 
  \emph{d}-dimensional multivariate normal with covariance matrix \eqn{\hat{A}}. 
  Typically, \eqn{\hat{A}} is taken as the difference between the sample covariance and 
  the covariance of fitted log-concave maximum likelihood density. Therefore, this estimator 
  matches both the empirical mean and empirical covariance.
  
  The estimate is evaluated numerically either by Gaussian quadrature in two dimensions, or in 
  higher dimensions, via a combinatorial method proposed by \cite{Grundmann and Moeller (1978)}.
  Details of the computational aspects can be found in \cite{Chen and Samworth (2011)}. In one 
  dimension, explicit expression can be derived. See 
  \code{\link[logcondens:logcondens-package]{logcondens}} for more information.

  For examples, see \code{\link{mlelcd}}
}

\value{ A \code{vector} of smoothed log-concave maximum likelihood estimate
   values, as evaluated at the points \code{x}.}

\author{
  Yining Chen

  Madeleine Cule
  
  Robert Gramacy
  
  Richard Samworth
}

\references{
 Chen, Y. and Samworth, R. J. (2013) 
 \emph{Smoothed log-concave maximum likelihood estimation with applications} 
 Statist. Sinica, 23, 1373-1398. \url{https://arxiv.org/abs/1102.1191v4}

 Grundmann, A. and Moeller, M. (1978)
 \emph{Invariant Integration Formulas for the N-Simplex by Combinatorial Methods} 
 SIAM Journal on Numerical Analysis, Volume 15, Number 2, 282-290.

}

\seealso{ \code{\link{dlcd}}, \code{\link{hatA}}, \code{\link{mlelcd}} }

\keyword{multivariate}
\keyword{nonparametric}
\keyword{smoothing}
