#Function to evaluate the given density at a point (x,y)


'is.in' <- function(po,norm,offset) {
  invec <- NULL
  for (i in 1:nrow(po)) {
  inside <- apply(norm * t(po[i,] - t(offset)),1,sum)
  invec <- c(invec,sum(inside > 0) == 0)
}
  return(invec)
}

'eval2' <- function(po, out,uselog=FALSE) {
  isvec <- is.in(po,out$outnorm, out$outoffset)
  est <- rep(-Inf,nrow(po))
  for (i in 1:nrow(po)) {
    if(isvec[i])    est[i] <- min(apply(out$b,1,'%*%',po[i,])-out$beta)
  }
  if(uselog) return(est)
  else return(exp(est))
}




'lcd.eval' <- function(out, po, uselog=FALSE) {
  if(class(out) != "LogConcDEAD") {
    stop("error: out must be of class LogConcDEAD")
  }
  d <- ncol(out$x)
  if(is.vector(po) && length(po)==d) {
    po <- matrix(po,ncol=d)
  }
  if(is.matrix(po) && ncol(po)==d)
  {
    return(eval2(po,out,uselog))
  }
  else stop("error: po must be a vector, or a numeric matrix with ",d," columns")
}

 

