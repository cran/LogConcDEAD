plot.LogConcDEAD <- function (x, uselog = FALSE, type = "ic", addp = TRUE, 
    gridlen = 100, g = NULL, marg = NULL, g.marg = NULL, ...) 
{
  d <- ncol(x$x)

  ## one dimensional data
  if (d == 1) {
    y <- x$x[, 1]
    o <- order(y)
    type <- "l"
    if (uselog) { 
      plot(y[o], x$logMLE[o], type = type, ylab = "Log density estimate", 
           xlab = "X", ...)
    } else {
      plot(y[o], exp(x$logMLE[o]), ylab = "Density estimate", 
              xlab = "X", type = type, ...)
    }
  }

  ## marginals
  else if (!is.null(marg)) {

    ##If we have a valid marginal
    if (is.element(marg, 1:d)) {

      ##Check whether already calculated and, if not, calculate
      if (is.null(g.marg)) {
        g.marg <- lcd.marg(x,marg=marg)
      }
      
      if(uselog) {
        plot(g.marg$xo, log(g.marg$marg), type = "l", xlab = paste("X",marg),ylab = "log estimated marginal density", main=paste("Marginal for X",marg),...) }
      else {
        plot(g.marg$xo, g.marg$marg, type = "l", xlab = paste("X",marg), 
                 ylab = "estimated marginal density", main=paste("Marginal for X",marg),...)
      }
    }    
    else stop(cat("Marginal should be one of 1, ...",d,"\n"))
  }
  else {

    if (d > 2) {
      stop("It is not possible to plot for d>2. Marginal estimates may be plotted by setting the parameter marg")
    }
      y <- x$x
      z <- x$logMLE
      if (is.null(g)) {
        if (!require("akima", quietly = TRUE)) 
          stop("You need to install the akima package")
        g <- lcd.interp(x, gridlen = gridlen)
      }
      
      if(uselog) {
         mytitle <- "Log density estimate"
      } else {
        g$z <- exp(g$z)
        g$z[is.na(g$z)] <- 0
        mytitle <- "Density estimate"
      }
      if (type == "p") 
        persp(g, zlab = mytitle, xlab = "X_1", 
                ylab = "X_2", ...)
      else if (type == "i") 
        image(g, col = terrain.colors(128), main=mytitle,...)
      else if (type == "c") 
        contour(g, main=mytitle,...)
      else if (type == "ic") {
        image(g, col = terrain.colors(128),main=mytitle, ...)
        contour(g, add = TRUE, ...)
      }
      else if (type == "r") {
        if (!require("rgl", quietly = TRUE)) 
          stop("you need to install the rgl package")
        r <- range(as.numeric(g$z[is.finite(g$z)]))
        if (uselog) 
          g$z <- 0.2 * (r[2] - r[1]) * g$z
        else g$z <- 100* (r[2] - r[1]) * g$z
        zlim <- range(g$z[!is.na(g$z)])
        zlen <- zlim[2] - zlim[1] + 1
        colorlut <- terrain.colors(zlen)
        col <- colorlut[g$z - zlim[1] + 1]
        open3d()
        if (!uselog) {
          surface3d(g$x, g$y, g$z, color = col, back = "lines")
          decorate3d(xlim = range(y[, 1]), ylim = range(y[, 
                                             2]), zlim = range(g$z), zlab = "Density", ...)
        }
        else {
          decorate3d(xlim = range(y[, 1]), ylim = range(y[, 
                                             2]), zlim = range(g$z[!is.na(g$z)]), zlab = "Log density", 
                     ...)
          rgl.surface(g$x, g$y, g$z, coords = c(1, 3, 2), 
                      color = col, back = "lines")
        }
        if (addp) {
          if (uselog) 
            plot3d(x$x[, 1], x$x[, 2], 0.2 * (r[2] - r[1]) * 
                   z, pch = 4, size = 2, add = TRUE, col = "black")
          else plot3d(x$x[, 1], x$x[, 2], 100 * (r[2] - r[1]) * 
                      exp(z), pch = 4, size = 2, add = TRUE, col = "black")
        }
      }
      else {
        stop("type should be one of r, p, i, c, or ic")
      }
      if (addp && (type == "i" || type == "c" || type == 
                   "ic")) 
        points(y, ...)
    }
  }
