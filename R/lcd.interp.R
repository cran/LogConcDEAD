#Function to do the 2d interpolation (on the log scale)
'lcd.interp' <-
  function(out, gridlen=100, span=0.5, ...)
{

  if(!require("akima",quietly=TRUE)) stop("you need to install the akima package")
x <- out$x[,1]
y <- out$x[,2]
z <- out$logMLE
  
xo <- seq(min(x), max(x), length=gridlen)
yo <- seq(min(y), max(y), length=gridlen)

chull <- out$chull
nsimplices <- nrow(chull)

g <- matrix(NA,nrow=gridlen,ncol=gridlen)
for (i in 1:nsimplices)
  {
  verts <- chull[i,]
  g <- pmax(interp(x=c(x[verts],mean(x[verts])),y=c(y[verts],mean(y[verts])),z=c(z[verts],mean(z[verts])),xo=xo,yo=yo)$z,g,na.rm=TRUE)
  }
return(list(x=xo,y=yo,z=g))
}
