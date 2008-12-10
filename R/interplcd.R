#Function to do the 2d interpolation (on the log scale)

'interplcd' <-
  function(lcd, gridlen=100, span=0.5, ...)
{
  if(!require("akima",quietly=TRUE)) stop("you need to install the akima package")
  x <- lcd$x[,1]
  y <- lcd$x[,2]
  z <- lcd$logMLE
  
  xo <- seq(min(x), max(x), length=gridlen)
  yo <- seq(min(y), max(y), length=gridlen)

  chull <- lcd$triang
  nsimplices <- nrow(chull)

  g <- matrix(NA,nrow=gridlen,ncol=gridlen)

  for (i in 1:nsimplices)
  {
    verts <- chull[i,]
    if( sum( !is.finite( z[verts] ) ) == 0 ) {
      g <-
  pmax(interp(x=c(x[verts],mean(x[verts])),y=c(y[verts],mean(y[verts])),z=c(z[verts],mean(z[verts])),xo=xo,yo=yo)$z,g,na.rm=TRUE)
    }
  }
  return(list(x=xo,y=yo,z=g))
}

'lcd.interp' <- function( lcd, gridlen=100, span=0.5, ... ) {
  warning("lcd.interp is deprecated.  Use interplcd instead")
  return( interplcd( lcd, gridlen, span, ... ) )
}
