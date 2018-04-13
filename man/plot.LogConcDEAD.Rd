\name{plot.LogConcDEAD}
\alias{plot.LogConcDEAD}

\title{Plot a log-concave maximum likelihood estimator} \description{
\code{\link[graphics:plot]{plot}} method for class \code{"LogConcDEAD"}.  Plots of various
types are available for 1- and 2-d data. For dimension greater than 1,
plots of axis-aligned marginal density estimates are available.}

\usage{
\method{plot}{LogConcDEAD}(x, uselog=FALSE, type="ic", addp=TRUE,
  drawlabels=TRUE, gridlen=400, g, marg, g.marg, main, xlab, ylab, ...)
}

\arguments{
  \item{x}{Object of class \code{"LogConcDEAD"} (typically output from \code{\link{mlelcd}})}
  \item{uselog}{Scalar \code{logical}: should the plot be on the log scale?}
  \item{type}{Plot type: \code{"p"} perspective, \code{"c"} contour, \code{"i"} image, \code{ic}
    image and contour, \code{r} using \pkg{\link[rgl:rgl-package]{rgl}} (the best!)}
  \item{addp}{Scalar \code{logical}: should the data points be plotted? (as black dots on the surface for \eqn{d \geq 2}{d >= 2}; as circles for \eqn{d=1})}
  \item{drawlabels}{Scalar \code{logical}: should labels be added to
  contour lines? (only relevant for \code{type}s \code{"ic"} and \code{"c"})}
  \item{gridlen}{Integer scalar indicating the number of points at which the maximum likelihood estimator is evaluated in each
    dimension}
  \item{g}{(optional) a \code{matrix} of density estimate values (the result of a call to  \code{\link{interplcd}}).  If many plots of a single dataset are
  required, it may be quicker to compute the grid using
  \code{\link{interplcd}(x)} and pass the result to \code{plot}}
  \item{marg}{If non-\code{NULL}, this scalar integer determines which marginal should be
  plotted (should be between \eqn{1} and \eqn{d})}
  \item{g.marg}{If \code{g} is non-\code{NULL}, can contain a
  \code{vector} of marginal density estimate values (the output of \code{\link{interpmarglcd}}).  If many plots of a single dataset
  are required, it may be quicker to compute the marginal values to
  compute marginal values using \code{\link{interpmarglcd}}
  and pass the result to \code{plot}}
\item{main}{Title}
\item{xlab}{x-axis label}
\item{ylab}{y-axis label}
\item{...}{Other arguments to be passed to the generic
  \code{\link[graphics:plot]{plot}} method}
}
\details{
  The density estimate is evaluated on a grid of points using the
  \code{\link{interplcd}} function. If several plots are required, this
  may be computed separately and passed to \code{plot} using the
  \code{g} argument.
  
  For two dimensional data, the default plot type is \code{"ic"},
  corresponding to \code{\link{image}} and \code{\link{contour}} plots.
  These may be obtained separately using plot type \code{"i"} or \code{"c"}
  respectively.  Where available, the use of plot type \code{"r"} is
  recommended.  This uses the \pkg{\link[rgl:rgl-package]{rgl}} package
  to produce a 3-d plot that may be rotated by the user.  The option
  \code{"p"} produces perspective plots.
  
  For data of dimension at least 2, axis-aligned marginals may be
  plotted by setting the \code{marg} argument.  This integrates the
  estimated density over the remaining dimensions. If several plots are
  required, the estimate may be computed using the function
  \code{\link{interpmarglcd}} and passed using the argument
  \code{g.marg}.

  Where relevant, the colors were obtained from the function
  \code{\link[colorspace]{heat_hcl}} in the package
  \pkg{colorspace}. Thanks to
  Achim Zeileis for this suggestion.

  For examples, see \code{\link{mlelcd}}.
}


\author{Madeleine Cule

  Robert B. Gramacy

  Richard Samworth
  
  Yining Chen}

\seealso{\code{\link{mlelcd}}, \code{\link{interplcd}}, \code{\link{interpmarglcd}}, \code{\link[colorspace]{heat_hcl}} }



\keyword{ dplot }
\keyword{ hplot }
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ iplot }
\keyword{ dynamic }
