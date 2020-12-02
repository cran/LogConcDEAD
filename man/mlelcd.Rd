\name{mlelcd}
\alias{mlelcd}

\title{Compute the maximum likelihood estimator of a log-concave density}

\description{
  Uses Shor's \eqn{r}{r}-algorithm to compute the maximum
 likelihood estimator of a log-concave density based on an
 i.i.d. sample.  The estimator is uniquely determined by its value at the data
 points. The output is an object of class \code{"LogConcDEAD"} which
  contains all the information needed to plot the estimator using the
  \code{\link[=plot.LogConcDEAD]{plot}} method, or to evaluate it using
  the function \code{\link{dlcd}}.
 }

 \usage{
mlelcd(x, w=rep(1/nrow(x),nrow(x)), y=initialy(x),
  verbose=-1, alpha=5, c=1, sigmatol=10^-8, integraltol=10^-4,
  ytol=10^-4, Jtol=0.001, chtol=10^-6)
 }

\arguments{
  \item{x}{Data in \eqn{R^d}, in the form of an \eqn{n
      \times d}{n x d} numeric \code{matrix}}
  \item{w}{Vector of weights \eqn{w_i} such that the computed estimator
  maximizes \deqn{\sum_{i=1}^n w_i \log f(x_i)}{w[1] log f(x[1,]) + ... +
  w[n] log f([x,n])}} subject to the restriction that \eqn{f} is
  log-concave. The default is  \eqn{\frac{1}{n}}{1/n} for all \eqn{i},
  which corresponds to i.i.d. observations.   
  \item{y}{Vector giving starting point for the \eqn{r}-algorithm.
    If none given, a kernel estimate is used.}
  \item{verbose}{
    \itemize{
      \item -1: (default) prints nothing
      \item 0: prints warning messages
      \item \eqn{n>0}: prints summary information every \eqn{n}
      iterations
    }
  }
  \item{alpha}{Scalar parameter for \pkg{SolvOpt}}
  \item{c}{Scalar giving starting step size}
  \item{sigmatol}{Real-valued scalar giving one of the stopping
  criteria: Relative change in \eqn{\sigma}{sigma} must be below
  \code{sigmatol} for algorithm to terminate. (See Details)}
  \item{ytol}{Real-valued scalar giving on of the stopping criteria: Relative change in \eqn{y} must be
  below \code{ytol} for algorithm to terminate. (See
    Details)}
  \item{integraltol}{Real-valued scalar giving one of the stopping
    criteria:  \eqn{| 1 - \exp(\bar{h}_y) |}{ |exp(h_y) - 1|} must be below
  \code{integraltol} for algorithm to terminate. (See Details)}
  \item{Jtol}{Parameter controlling when Taylor expansion is used in
    computing the function \eqn{\sigma}{sigma}}
  \item{chtol}{Parameter controlling convex hull computations}
}

\details{
  The log-concave maximum likelihood density estimator based on data
  \eqn{X_1, \ldots, X_n}{X_1, ..., X_n} is the function that maximizes
  \deqn{\sum_{i=1}^n w_i \log f(X_i)}{(w_1 log f(X_1) + ... +
    w_n log f(X_n))} subject to the constraint that \eqn{f} is
  log-concave. For i.i.d.~data, the weights \eqn{w_i} should be
  \eqn{\frac{1}{n}}{1/n} for each \eqn{i}.
  
  This is a function of the form \eqn{\bar{h}_y}{h_y} for some \eqn{y
    \in R^n}{y in R^n}, where \deqn{\bar{h}_y(x)
    = \inf \lbrace h(x) \colon h \textrm{ concave }, h(x_i) \geq y_i
    \textrm{ for } i = 1, \ldots, n \rbrace.}{h_y(x) =
    inf{h(x): h concave, h(x_i) >= y_i for i = 1, ..., n}.}

  Functions of this form may equivalently be specified by dividing
  \eqn{C_n}, the convex hull of the data, into simplices \eqn{C_j} for
  \eqn{j \in J} (triangles in 2d, tetrahedra in 3d etc), and setting
  \deqn{f(x) = \exp\{b_j^T x - \beta_j\}}{f(x) = exp{b_j^T x - beta_j}}
  for \eqn{x \in C_j}{x in C_j}, and \eqn{f(x) = 0} for \eqn{x \notin
  C_n}{x not in C_n}.
  
  This function uses Shor's \eqn{r}-algorithm (an iterative
  subgradient-based procedure) to minimize over vectors \eqn{y}{y} in
  \eqn{R^n} the function \deqn{\sigma(y) = -\frac{1}{n} \sum_{i=1}^n y_i
  + \int \exp(\bar{h}_y(x)) \, dx.}{sigma(y) = -1/n (y_1 + ... + y_n) +
  int (h_y(x)) dx.}  This is equivalent to finding the log-concave
  maximum likelihood estimator, as demonstrated in \cite{Cule, Samworth
  and Stewart (2008)}.
  
  An implementation of Shor's \eqn{r}-algorithm based on \pkg{SolvOpt}
is used.

  Computing \eqn{\sigma}{sigma} makes use of the \pkg{qhull} library.
Code from this C-based library is copied here as it is not currently 
possible to use compiled code from another library.  
For points not in general position, this requires a Taylor expansion of 
\eqn{\sigma}{sigma}, discussed in \cite{Cule and D\"umbgen (2008)}.
}

\value{
  An object of class \code{"LogConcDEAD"}, with the following
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

\references{
Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T. (1996)
\emph{The Quickhull algorithm for convex hulls}
  ACM Trans. on Mathematical Software, 22(4) p.469-483
\url{http://www.qhull.org}

Cule, M. L. and D\"umbgen, L. (2008) On an auxiliary function for
log-density estimation, University of Bern technical report. \url{https://arxiv.org/abs/0807.4719}

Cule, M. L., Samworth, R. J., and Stewart, M. I. (2010)
\emph{Maximum likelihood estimation of a log-concave density},
  Journal of the Royal Statistical Society, Series B, 72(5) p.545-607.

Kappel, F. and Kuntsevich, A. V. (2000)
\emph{An implementation of Shor's r-algorithm},
  Computational Optimization and Applications, Volume 15, Issue 2, 193-205.

Shor, N. Z. (1985)
\emph{Minimization methods for nondifferentiable functions},
  Springer-Verlag
 }

\note{
  For one-dimensional data, the active set algorithm of
  \code{\link[logcondens:logcondens-package]{logcondens}} is faster, and may be
  preferred.

  The authors gratefully acknowledge the assistance of Lutz Duembgen
  at the University of Bern for his insight into the objective function
  \eqn{\sigma}{sigma}.

  Further references, including definitions and background material, may be found in \cite{Cule, Samworth and Stewart (2010)}.
}

\author{
  Madeleine Cule 
  
  Robert B. Gramacy

  Richard Samworth
  
  Yining Chen
}

\examples{
## Some simple normal data, and a few plots

x <- matrix(rnorm(200),ncol=2)
lcd <- mlelcd(x)
g <- interplcd(lcd)
par(mfrow=c(2,2), ask=TRUE)
plot(lcd, g=g, type="c")
plot(lcd, g=g, type="c", uselog=TRUE)
plot(lcd, g=g, type="i")
plot(lcd, g=g, type="i", uselog=TRUE)

## 2D interactive plot (need rgl package, not run here)
#  plot(lcd, type="r")


## Some plots of marginal estimates
par(mfrow=c(1,1))
g.marg1 <- interpmarglcd(lcd, marg=1)
g.marg2 <- interpmarglcd(lcd, marg=2)
plot(lcd, marg=1, g.marg=g.marg1)
plot(lcd, marg=2, g.marg=g.marg2) 

## generate some points from the fitted density
## via independent rejection sampling
generated1 <- rlcd(100, lcd)
colMeans(generated1)
## via Metropolis-Hastings algorithm
generated2 <- rlcd(100, lcd, "MH")
colMeans(generated2)

## evaluate the fitted density
mypoint <- c(0, 0)
dlcd(mypoint, lcd, uselog=FALSE)
mypoint <- c(1, 0)
dlcd(mypoint, lcd, uselog=FALSE)

## evaluate the marginal density
dmarglcd(0, lcd, marg=1)
dmarglcd(1, lcd, marg=2)

## evaluate the covariance matrix of the fitted density
covariance <- cov.LogConcDEAD(lcd)

## find the hat matrix for the smoothed log-concave that
## matches empirical mean and covariance
A <- hatA(lcd)

## evaluate the fitted smoothed log-concave density
mypoint <- c(0, 0)
dslcd(mypoint, lcd, A)
mypoint <- c(1, 0)
dslcd(mypoint, lcd, A)

## generate some points from the fitted smoothed log-concave density
generated <- rslcd(100, lcd, A)
}

\seealso{
  \code{\link[logcondens:logcondens-package]{logcondens}}, 
  \code{\link{interplcd}},  \code{\link{plot.LogConcDEAD}}, 
  \code{\link{interpmarglcd}}, \code{\link{rlcd}}, \code{\link{dlcd}},

  \code{\link{dmarglcd}}, \code{\link{cov.LogConcDEAD}}
}
\keyword{ nonparametric }
\keyword{ multivariate }
