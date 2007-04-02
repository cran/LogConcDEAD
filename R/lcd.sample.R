'mle.samples' <- function(prob,z,A,alpha,y1,dim,nsample=1) {
  sample <- NULL
  nrows <- length(prob)
  simp <- sample(1:nrows,nsample,prob=prob,replace=TRUE)
  for (i in 1:nsample) {
    x <- NULL
    while(is.null(x)) {
    w <- sort(runif(dim))
    w <- c(w[1],diff(w))
    fw <- exp(z[simp[i],]%*%w+y1[simp[i]])
    maxz <- exp(y1[simp[i]] + max(c(0,z[simp[i],])))
    u <- runif(1)
    if(u<fw/maxz) x <-(A[[simp[i]]])%*%w + alpha[simp[i],]
   # else print("no luck\n")

  }
    sample <- cbind(sample,x)
    }
  return(t(sample)) }

'lcd.sample' <- function(out,nsample=1) {
  chull <- out$chull
  x <- out$x
  y <- out$logMLE
  nrows <- nrow(chull)
  d <- ncol(x)
  A <- vector(mode="list",length=nrows)
  alpha <- matrix(0,nrow=nrows,ncol=d)
  y1 <- vector(length=nrows)
  prob <- vector(length=nrows)
  z <- matrix(0,nrow=nrows,ncol=d)
  for (j in 1:nrows) {
    vertices <- x[chull[j,],]
    A[[j]] <- t(vertices[-1,] - rep(1,d) %o% vertices[1,])
    z[j,] <- y[chull[j,-1]]-y[chull[j,1]]
    alpha[j,] <- vertices[1,]
    y1[j] <- y[chull[j,1]]
    zProd <- rep(0,d)
    for(r in 1:d) zProd[r] <- prod(z[j,r]-z[j,-r])
    prob[j] <- abs(det(A[[j]]))*exp(y1[j])*sum((exp(z[j,])-1)/(z[j,]*zProd))
  }
  prob <- prob/sum(prob)
mle.samples(prob,z,A,alpha,y1,dim=d,nsample)

}
