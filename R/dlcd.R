#Function to evaluate the given density at a point (x,y)


'is.in' <- function(po,norm,offset) {
  invec <- NULL
  for (i in 1:nrow(po)) {
  inside <- apply(norm * t(po[i,] - t(offset)),1,sum)
  invec <- c(invec,sum(inside > 0) == 0)
}
  return(invec)
}

'eval2' <- function(po, lcd,uselog=FALSE) {
  isvec <- is.in(po,lcd$outnorm, lcd$outoffset)
  est <- rep(-Inf,nrow(po))
  for (i in 1:nrow(po)) {
    if(isvec[i])    est[i] <- min(apply(lcd$b,1,'%*%',po[i,])-lcd$beta)
  }
  if(uselog) return(est)
  else return(exp(est))
}




'dlcd' <- function(x, lcd, uselog=FALSE) {
  if(class(lcd) != "LogConcDEAD") {
    stop("error: lcd must be of class LogConcDEAD")
  }
  d <- ncol(lcd$x)
  if(is.vector(x) && length(x)==d) {
    x <- matrix(x,ncol=d)
  }
  if(is.matrix(x) && ncol(x)==d)
  {
    return(eval2(x,lcd,uselog))
  }
  else stop("error: x must be a vector, or a numeric matrix with ",d," columns")
}

###Retained for compatibility
'lcd.eval' <- function(lcd, po, uselog=FALSE) {
  warning("lcd.eval is deprecated.  Use dlcd instead")
  return(dlcd( po, lcd, uselog ) )
}
 

