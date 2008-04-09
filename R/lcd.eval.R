#Function to evaluate the given density at a point (x,y)

'eval.vect' <- function(yvect, x, b, beta, uselog=FALSE) {
  y <- as.matrix(yvect)
  return(as.vector(apply(y,1,eval.est,x=x,b=b,beta=beta,uselog=uselog)))}

'eval.est' <- function(y, x, b, beta, uselog=FALSE) {
  point <- c(x,y)
  logdens <- (min(apply(b,1,'%*%',point)-beta))
  if (uselog) return (logdens)
  else return (exp(logdens))
}

'eval.est2'  <- function(x, y, b, beta, uselog=FALSE) {
  point <- c(x,y)
  logdens <- (min(apply(b,1,'%*%',point)-beta))
  if (uselog) return (logdens)
  else return (exp(logdens)) }


'eval.vect2' <- function(xvect, y, b, beta, uselog=FALSE) {
  x <- as.matrix(xvect)
  return(as.vector(apply(x,1,eval.est,y=y,b=b,beta=beta,uselog=uselog)))}


#Find the limits (at the moment, just in terms of X_2 for a particular X_1)
'find.limits' <- function(x, X, chull) {
  nrows <- nrow(chull)
  lim <- NULL

  for (i in 1:nrows) {
    y <- X[chull[i,1],]
    z <- X[chull[i,2],]
    if((x >= min(y[1],z[1]) && x <= max(y[1],z[1]))) {
    lim <- c(lim,y[2] + (x-y[1])*(z[2]-y[2])/(z[1]-y[1]))
    }
  }
  if(is.null(lim)) lim <- c(0,0)
  if(length(lim)==1) lim <- c(lim,lim)
  return(sort(lim))
}


'find.limits2' <- function(x, X,chull) {
  nrows <- nrow(chull)
  lim <- NULL

  for (i in 1:nrows) {
    y <- X[chull[i,1],]
    z <- X[chull[i,2],]
    tmp <- y[1] + (x-y[2])*(z[1]-y[1])/(z[2]-y[2])
    if((tmp >= min(y[1],z[1])) && (tmp <= max(y[1],z[1]))) {
      lim <- c(lim,tmp)
    }
  }
  if(is.null(lim)) lim <- c(0,0)
  if(length(lim)==1) lim <- c(lim,lim)
  return(sort(lim))
}



#Evaluate the marginal density at a point (currently only does X_1...)
'eval.marg' <- function(x,out) {
  ##if(!require("geometry",quietly=TRUE)) stop("you need to install the geometry package")
  beta <- out$beta
  b <- out$b
  X <- out$x
  chull <- convhullnew(X)
  limits <- find.limits(x,X,chull)
  return(integrate(eval.vect,lower=limits[1],upper=limits[2],x=x,b=b,beta=beta)$value)
}

'eval.marg2' <- function(y,out) {
  ##if(!require("geometry",quietly=TRUE)) stop("you need to install the geometry package")
  beta <- out$beta
  b <- out$b
  X <- out$x
  chull <- convhullnew(X)
  limits <- find.limits2(y,X,chull)
  return(integrate(eval.vect2,lower=limits[1],upper=limits[2],y=y,b=b,beta=beta)$value)
} 

## Finds the marginal density of X_1
'lcd.marg'<- function(out, gridlen=100) {
  X <- out$x
  xo <- matrix(seq(min(X[,1]),max(X[,1]),length=gridlen))
  marg = apply(xo,1,eval.marg,out=out)
  return(list(marg=marg,xo=xo))
}

  
'lcd.marg2' <- function(out,gridlen=100) {
  X <- out$x
  yo <- matrix(seq(min(X[,2]),max(X[,2]),length=gridlen))
  marg = apply(yo,1,eval.marg2,out=out)
  return(list(marg=marg,yo=yo)) }

'lcd.eval' <- function(points ,out, uselog=FALSE) {
  n <- nrow(out$x)
  d <- ncol(out$x)

  if(is.vector(points) && length(points)==d) {
    points <- t(as.matrix(points)) }
  ## check the dimensions of points
  if (is.matrix(points) && ncol(points)==d) {
  return (apply(points,1,lcd.eval.point,out=out,uselog=uselog))
}
  else stop("error: points must be a vector, or a numeric matrix with ",d," columns")
}

'lcd.eval.point' <- function(point, out, uselog=FALSE) {
  chull <- out$chull
  verts <- out$verts
  vertsoffset <- out$vertsoffset
  b <- out$b 
  beta <- out$beta
  d <- ncol(out$x)
  lambda <- t(apply(verts,1,'%*%',point)) - vertsoffset
  which <- ((apply((lambda>=0)&(lambda<=1),1,sum)==d)*(apply(lambda,1,sum)<=1)==1)
  if(sum(which)==0)  {
    if(uselog) return(-Inf)
    else return(0) }
  else {
    which <- order(which,decreasing=TRUE)[1]
    if(uselog)
      return(b[which,]%*%point - beta[which])
    else
      return(exp(b[which,]%*%point - beta[which]))
  }
}
