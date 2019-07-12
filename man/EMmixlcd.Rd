\name{EMmixlcd}
\alias{EMmixlcd}

\title{Estimate the mixture proportions and component densities using EM algorithm}

\description{
  Uses EM algorithm to estimate the mixture proportions and the component
  densities. The output is an object of class \code{"lcdmix"} which
  contains mixture proportions at each observation and all the information of the 
  estimated component densities.
}

\usage{
  EMmixlcd( x, k = 2, y, props, epsratio=10^-6, max.iter=50,
            epstheta=10^-8, verbose=-1 )
}

\arguments{
  \item{x}{Data in \eqn{R^d}, in the form of an \eqn{n
      \times d}{n x d} numeric \code{matrix}}
  \item{k}{The number of components, equals 2 by default}   
  \item{y}{An \eqn{n \times k}{n x k} numeric \code{matrix} giving the starting values for the EM 
      algorithm. If none given, a hierachical Gaussian clustering model is used. To reduce the 
      computational burden while allowing sufficient flexibility for the EM algorithm, it is recommended 
      to leave this argument unspecified.}
  \item{props}{Vector of length \eqn{k} containing the starting value of proportions.
      If none given, a hierachical Gaussian clustering model is used. To reduce the 
      computational burden while allowing sufficient flexibility for the EM algorithm, 
      it is recommended to leave this argument unspecified.}
  \item{epsratio}{EM algorithm will terminate if the increase in the proportion of the likelihood is 
      less than this specified ratio. Default value is \eqn{10^{-6}}. }
  \item{max.iter}{The maximum number of iterations for the EM algorithm}
  \item{epstheta}{\eqn{epstheta/n} is the thresold of the weight below which data point is 
      discarded from the cluster. This quantity is introduced to increase the computational efficiency
      and stability.}
  \item{verbose}{
    \itemize{
      \item -1: (default) prints nothing
      \item 0: prints warning messages
      \item \eqn{>0}: prints summary information every \eqn{n} iterations
    }
  }
}

\details{
  An introduction to the Em algorithm can be found in \cite{McLachlan and Krishnan (1997)}.
  Briefly, given the current estimates of the mixture proportions and component densities, we first 
  update the estimates of the mixture prroportions. We then update the estimates of the component 
  densities by using \code{\link{mlelcd}}. In fact, the incorporation of the weights in the 
  maximization process in \code{\link{mlelcd}} presents no additional complication.

  In our case, because of the computational intensity of the method, we first cluster the points 
  according to ta hierarchical Gaussian clustering model and then iterate the EM algorithm until the 
  increase in the proportion of the likelihood is less than a pre-specified quantity at each step.

  More technical details can be found in \cite{Cule, Samworth and Stewart(2010)}
}

\value{
  An object of class \code{"lcdmix"}, with the following components:
  \item{x}{Data copied from input (may be reordered)}
  \item{logf}{An \eqn{n \times k}{n x k} \code{maxtrix} of the log of the maximum likelihood estimate, 
    evaluated at the observation points for each component.}
  \item{props}{Vector containing the estimated proportions of components}
  \item{niter}{Number of iterations of the EM algorithm}
  \item{lcdloglik}{The log-likelihood after the final iteration}
}

\references{

Cule, M. L., Samworth, R. J., and Stewart, M. I. (2010)
\emph{Maximum likelihood estimation of a log-concave density},
  Journal of the Royal Statistical Society, Series B, 72(5) p.545-607.

McLachlan, G. J. and Krishnan, T. (1997) 
\emph{The EM Algorithm and Extensions}, New York: Wiley.

}


\author{
  Yining Chen
  
  Madeleine Cule 
  
  Robert B. Gramacy

  Richard Samworth}

\examples{
##Simple bivariate normal data
  set.seed( 1 )
  n = 15
  d = 2
  props=c( 0.6, 0.4 )
  shift=2
  x <- matrix( rnorm( n*d ), ncol = d )
  shiftvec <- ifelse( runif( n ) > props[ 1 ], 0, shift )
  x[,1] <- x[,1] + shiftvec
  EMmixlcd( x, k = 2, max.iter = 2)
}

\seealso{
  \code{\link[mclust:hc]{mclust}},
  \code{\link[logcondens:logcondens-package]{logcondens}},
  \code{\link{plot.LogConcDEAD}},\code{\link{mlelcd}}, \code{\link{dlcd}}
}
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ EM }
\keyword{ classification }
