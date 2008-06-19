##This function is internal and computes the integral of
## exp(b^T x - beta) over
## T \cup L (w.r.t. Lebesgue measure on an appropriate subspace)
## where T is the simplex with vertices given by verts
## and L is the subspace X[,keep] = lambda

'integratemarg' <- function(keep, lambda, b, beta, verts, d=2, eps=10^-6) {
  vals <- verts[,keep] - lambda
  left <- vals < -eps
  right <- vals > eps
  on <- abs(vals) <= eps
  ans <- 0
  new <- NULL
  new <- rbind(new,verts[on,])
  for (i in which(left)) {
    for (j in which(right)) {
      new <- rbind(new, verts[i,] + (lambda - verts[i,keep])/(verts[i,keep] - verts[j,keep]) *(verts[i,] - verts[j,]))
    }
  }

  if(nrow(new)==d) {
    tmp <- new[,-keep,drop=FALSE]
    y <- new%*%b - beta
    ans <- ans + J00AD(y,d-1)%*%abs(det(tmp[-1,,drop=FALSE]-rep(1,d-1)%*%tmp[1,,drop=FALSE]))
  }
  else if(nrow(new)>d) {
    tmp <- new[,-keep,drop=FALSE]
    trig <- matrix(delaunayn(tmp),ncol=d)
    for (i in 1:nrow(trig)) {
      y <- new[trig[i,],]%*%b - beta
      ans <- ans + J00AD(y,d-1)%*%abs(det(tmp[trig[i,-1],,drop=FALSE]-rep(1,d-1)%*%tmp[trig[i,1],,drop=FALSE]))
    }
  }
  return(ans)
}


##Computes one-dimensional marginals!

'lcd.marg' <- function(out, marg=1, gridlen=100) {
  if(class(out) != "LogConcDEAD") {
    stop("error: out must be of class LogConcDEAD")
  }
  x <- out$x
  d <- ncol(x)
  if(!is.element(marg,1:d)) {
    stop("error: marg must be one of 1 ... ",d)
  }
  triang <- out$triang
  xo <- seq(min(x[,marg]),max(x[,marg]),len=gridlen)
  val <- rep(0,gridlen)
  for (j in 1:nrow(triang)) {
    for (i in which(xo > min(x[triang[j,],marg]) & xo < max(x[triang[j,],marg]))) {
      val[i] <- val[i] + integratemarg(marg, xo[i], out$b[j,], out$beta[j], out$x[out$triang[j,],],d=d)
    }
  }
  return(list(xo=xo,marg=val))
}

'lcd.marg.eval' <- function(out, point=0, marg=1) {
  if(class(out) != "LogConcDEAD") {
    stop("error: out must be of class LogConcDEAD")
  }
  x <- out$x
  d <- ncol(x)
  if(!is.element(marg, 1:d)) {
    stop("error: marg must be one of 1 ... ",d)
  }
  triang <- out$triang
  ans <- 0
  for (j in 1:nrow(triang)) {
    ans <- ans + integratemarg(marg, point
                               , out$b[j,], out$beta[j], out$x[out$triang[j,],],d=d)
  }
  return(ans)
}


##Some R code kindly provided by Lutz Duembgen to compute the integrals
'J00AD_appr' <- function(y,d){
  ## Approximate
  ## J(y)  =  exp(mean(y)) * J(y - mean(y))
  ## via Taylor expansion to third order of J(y - mean(y)).
  ybar <- mean(y)
  z <- y - ybar
  f0 <- prod(1:d)
  f2 <- f0*(d+1)*(d+2)
  f3 <- f2*(d+3)
  z2 <- sum(z^2)/2
  z3 <- sum(z^3)/3
  return(exp(ybar)*(1/f0 + z2/f2 + z3/f3))
}

'J00AD_ord' <- function(y,d,eps=0.001){
  ## Recursive implementation J(y) for ordered vector y of length d+1.
  ## If y[d+1] - y[1] < eps, a certain Taylor approximation is used.
  if(y[d+1] - y[1] < eps){
    return(J00AD_appr(y,d))		
  }
  else{
    if(d == 1){
      return((exp(y[1]) - exp(y[2]))/(y[1] - y[2]))
    }
    else{
      e1 <- J00AD_ord(y[1:d]    ,d-1, eps)
      e2 <- J00AD_ord(y[2:(d+1)],d-1, eps)
      return((e1 - e2)/(y[1] - y[d+1]))
    }
  }
}

'J00AD' <- function(y,d,eps=0.001){
  ## Computes the function J(y) for a vector y of length d+1.
  y <- sort(y)
  return(J00AD_ord(y,d,eps))
}
