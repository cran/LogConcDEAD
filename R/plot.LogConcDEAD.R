plot.LogConcDEAD <- function (x, uselog = FALSE, method = "akima", itype = "p", addp = TRUE, 
    gridlen = 100, g = NULL, marg = NULL, g.marg = NULL, ...) 
{
    if (ncol(x$x) == 1) {
        y <- x$x[, 1]
        type <- "l"
        o <- order(y)
        if (uselog) 
            plot(y[o], x$logMLE[o], type = type, ylab = "Log density estimate", 
                xlab = "X", ...)
        else plot(y[o], exp(x$logMLE[o]), ylab = "Density estimate", 
            xlab = "X", type = type, ...)
    }
    else if (!is.null(marg) && ncol(x$x) == 2) {
        if (marg == 1) {
            if (is.null(g.marg)) 
                g.marg <- lcd.marg(x)
            plot(g.marg$xo, g.marg$marg, type = "l", xlab = "X_1", 
                ylab = "estimated marginal density", ...)
        }
        else if (marg == 2) {
            if (is.null(g.marg)) 
                g.marg <- lcd.marg2(x)
            plot(g.marg$yo, g.marg$marg, type = "l", xlab = "X_2", 
                ylab = "estimated marginal density", ...)
        }
        else stop("Marginal should be one of 1, or 2, and data should be 2d")
    }
    else {
        if (ncol(x$x) > 2) 
            stop("not done yet")
        y <- x$x
        z <- x$logMLE
        if (is.null(g)) {
            if (!require("akima", quietly = TRUE)) 
                stop("you need to install the akima package")
            g <- lcd.interp(x, gridlen = gridlen)
        }
        else g <- g
        if (!uselog) {
            g$z <- exp(g$z)
            g$z[is.na(g$z)] <- 0
        }
        if (itype == "p") 
            persp(g, zlab = "density estimate", xlab = "X_1", 
                ylab = "X_2", ...)
        else if (itype == "i") 
            image(g, col = terrain.colors(128), ...)
        else if (itype == "c") 
            contour(g, ...)
        else if (itype == "ic") {
            image(g, col = terrain.colors(128), ...)
            contour(g, add = TRUE, ...)
        }
        else if (itype == "r") {
            if (!require("rgl", quietly = TRUE)) 
                stop("you need to install the rgl package")
            r <- range(as.numeric(g$z[is.finite(g$z)]))
            if (uselog) 
                z2 <- 0.5 * (r[2] - r[1]) * g$z
            else z2 <- 500* (r[2] - r[1]) * g$z
            zlim <- range(z2[!is.na(z2)])
            zlen <- zlim[2] - zlim[1] + 1
            colorlut <- terrain.colors(zlen)
            col <- colorlut[z2 - zlim[1] + 1]
            open3d()
            if (!uselog) {
                surface3d(g$x, g$y, z2, color = col, back = "lines")
                decorate3d(xlim = range(y[, 1]), ylim = range(y[, 
                  2]), zlim = range(g$z), zlab = "Density", ...)
            }
            else {
                decorate3d(xlim = range(y[, 1]), ylim = range(y[, 
                  2]), zlim = range(g$z[!is.na(g$z)]), zlab = "Log density", 
                  ...)
                rgl.surface(g$x, g$y, z2, coords = c(1, 3, 2), 
                  color = col, back = "lines")
            }
            if (addp) {
                if (uselog) 
                  plot3d(x$x[, 1], x$x[, 2], 0.5 * (r[2] - r[1]) * 
                    z, pch = 4, size = 2, add = TRUE, col = "black")
                else plot3d(x$x[, 1], x$x[, 2], 500 * (r[2] - r[1]) * 
                  exp(z), pch = 4, size = 2, add = TRUE, col = "black")
            }
        }
        else {
            stop("itype should be one of r, p, i, c, or ic")
        }
        if (addp && (itype == "i" || itype == "c" || itype == 
            "ic")) 
            points(y, ...)
    }
}
