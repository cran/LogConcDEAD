\name{lcd.mle}
\alias{lcd.mle}

\title{Compute the maximum likelihood estimator of a log-concave density}
\description{
 This function uses Shor's \eqn{r}{r}-algorithm to compute the maximum
 likelihood estimator of a log-concave density based on an
 i.i.d. sample.  The estimator is determined by its value at data
 points.  }

\usage{lcd.mle (X, weights=rep(1/nrow(X),nrow(X)),
y=initialy(X), verbose=-1, alpha=5, c=1, sigmatol=10^-8,
integraltol=10^-4, ytol=10^-4, stepscale=5.1, stepscale2=2,
stepscale3=1.5, stepscale4=1.05, desiredsize=3.3, Jtol=0.001) }

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Data in \eqn{R^d}, in the form of an \eqn{n
      \times d}{n x d} numeric matrix}
  \item{weights}{weights \eqn{w_i} such that the computed estimator
  maximizes \eqn{\sum_{i=1}^n w_i f(X_i)}{w_1 f(X_1) + ... + w_n f(X_n)}}. The default is  \eqn{\frac{1}{n}}{1/n} for all \eqn{i}.   
  \item{y}{Starting point for the \eqn{r}-algorithm.  If none given, a kernel estimate is used.}
  \item{verbose}{
    \tabular{ll}{
      -1 \tab (default) prints nothing\cr 
      0  \tab  prints warning messages\cr
      \eqn{n>0} \tab  prints summary information every \eqn{n}
      iterations}
  }
  \item{alpha}{Scaling parameter for SolvOpt}
  \item{c}{Starting step size}
  \item{sigmatol}{Stopping criterion: relative change in
    \eqn{\sigma}{sigma} (see Details for details)}
  \item{ytol}{Stopping criterion: relative change in \eqn{y} (see
    Details for details)}
  \item{integraltol}{Stopping criterion: difference from 1 of integral
    of \eqn{\exp(\bar{h}_y)}{exp(h_y)} (see Details for details)}
  \item{stepscale,stepscale2,stepscale3,stepscale4}{Parameters for
    SolvOpt}
  \item{desiredsize}{Another parameter for SolvOpt (changing this is not
    recommended)}
  \item{Jtol}{Parameter controlling when Taylor expansion is used in
    computing the function \eqn{\sigma}{sigma}}
  }
\details{
  The MLE of a log-concave density is of the form
  \deqn{\bar{h}_y(x) = \inf \lbrace h(x) \colon h \textrm{ concave },
    h(X_i) \geq y_i \textrm{ for } i = 1, \ldots, n \rbrace}{h y(x) = inf{h(x): h concave, h(Xi) >= yi for i = 1, ..., n}}
  for some \eqn{y \in R^n}{y in R^n}.

  Functions of this form are piecewise linear on the convex hull of
  \eqn{X_1, \ldots, X_n}{X1, ..., Xn}, and of the form
  \deqn{\min \lbrace \langle b_j, x \rangle - \beta_j \rbrace}{min {<b_j, x> - beta_j}}
    for some vectors \eqn{b_j} and constants \eqn{\beta_j}{beta_j}.
  
  This function uses Shor's \eqn{r}-algorithm (an iterative
  subgradient-based procedure) to minimise over
  vectors \eqn{y}{y} in \eqn{R^n} the function
  \deqn{\sigma(y) = -\frac{1}{n} \sum_{i=1}^n y_i + \int
    \exp(\bar{h}_y(x)) \, dx.}{sigma(y) = -1/n (y_1 + ... + y_n) + int
    (h_y(x)) dx.}
  This is equivalent to finding the MLE.
  
  An implementation of Shor's \eqn{r}-algorithm based on
\emph{SolvOpt} is used.

  Computing \eqn{\sigma}{sigma} makes use of the \emph{qhull}
  (\url{http:\\www.qhull.org}) library, adapted from the \R
  implementation in \pkg{geometry}.  Code from this package is copied
  here as it is not currently possible to use compiled code from another
  package. 
}
\value{
  \item{x}{Data (reordered)}
  \item{w}{Weights}
  \item{logMLE}{Vector of
    the log of the MLE at the observation points}
  \item{NumberOfEvaluations}{Number of steps, number of function
    evaluations, number of subgradient evaluations.  If the SolvOpt
    algorithm fails, the first component will be an error message.}
  \item{MinSigma}{Minimum value of the objective function}
  \item{b}{see Details}
  \item{beta}{see Details}
  \item{chull}{Vector containing final triangulation of convex hull of data}
  \item{verts}{Details of triangulation for use in \code{\link{lcd.eval}}}
  \item{vertsoffset}{Details of triangulation for use in \code{\link{lcd.eval}}}
}
\references{
Cule, M. L., Samworth, R. J., and Stewart, M. I. (2007)
\emph{Computing the nonparametric maximum likelihood estimator of a
  log-concave density}
In preparation

Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T. (1996)
\emph{The Quickhull algorithm for convex hulls}
  ACM Trans. on Mathematical Software, 22(4) p.469-483
\url{http://www.qhull.org}

Kappel, F. and Kuntsevich, A. V. (2000)
\emph{An implementation of Shor's r-algorithm} Computational Optimization and Applications 15(2)
p.193-205

\url{http://www.uni-graz.at/imawww/kuntsevich/solvopt/}

N. Z. Shor (1985)
 \emph{Minimization methods for nondifferentiable functions}
 Springer-Verlag
}

\note{The authors gratefully acknowledge the assistance of Lutz Duembgen
  at the University of Bern for his insight into the objective function
  in \link{lcd.mle}.

Currently data not in general position are not supported.}

\author{Madeleine Cule \email{mlc40@cam.ac.uk}
  
  Robert B. Gramacy

  Richard Samworth}

\examples{
## some simple normal data
set.seed(101)
x <- matrix(rnorm(200), ncol=2)
out <- lcd.mle(x)
}
\seealso{
  \code{\link[<geometry>:<geometry>-package]{geometry}}
}
\keyword{ nonparametric }
\keyword{ multivariate }
