\name{interactive2D}
\alias{interactive2D}

\title{A GUI for classification in two dimensions using smoothed log-concave}

\description{
  Uses \code{\link[tkrplot:TkRplot]{tkrplot}} to create a GUI for two-class classification 
  in two dimensions using the smoothed log-concave maximum likelihood estimates
}

\usage{
  interactive2D(data, cl)
}

\arguments{
  \item{data}{Data in \eqn{R^2}, in the form of an \eqn{n
      \times 2}{n x 2} numeric \code{matrix}}
  \item{cl}{factor of true classifications of the data set}
}

\details{
  This function uses \code{\link[tkrplot:TkRplot]{tkrplot}} to create a GUI for two-class classification 
  in two dimensions using the smoothed log-concave maximum likelihood estimates. The construction of the 
  classifier is standard, and can be found in \cite{Chen and Samworth (2013)}. The slider controls the risk 
  ratio of two classes (equals one by default), which provides a way of demonstrating how the decision boundaries 
  change as the ratio varies. Observations from different classes are plotted in red and green respectively.
}

\value{
  A GUI with a slider
}

\references{

Chen, Y. and Samworth, R. J. (2013) 
 \emph{Smoothed log-concave maximum likelihood estimation with applications} 
 Statist. Sinica, 23, 1373-1398. \url{https://arxiv.org/abs/1102.1191v4}

Cule, M. L., Samworth, R. J., and Stewart, M. I. (2010)
\emph{Maximum likelihood estimation of a log-concave density},
  Journal of the Royal Statistical Society, Series B, 72(5) p.545-607.
}


\author{
  Yining Chen
  
  Madeleine Cule 
  
  Robert B. Gramacy

  Richard Samworth}

\examples{
## Simple bivariate normal data
## only works interactively, not run as a test example here
#  set.seed( 1 )
#  n = 15
#  d = 2
#  props=c( 0.6, 0.4 )
#  x <- matrix( rnorm( n*d ), ncol = d )
#  shiftvec <- ifelse( runif( n ) > props[ 1 ], 0, 1)
#  x[,1] <- x[,1] + shiftvec
#  interactive2D( x, shiftvec )
}

\seealso{
  \code{\link{dslcd}},\code{\link{mlelcd}}
}
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ classification }
