'lcd.mle' <-  function(X, initial=NULL, verbose=-1, alpha=5, c=1, sigmatol=10^-8, integraltol=10^-4, ytol=10^-4, stepscale=5.1, stepscale2=2,stepscale3=1.5,stepscale4=1.05,desiredsize=3.3) {

  
if (is.matrix(X)==FALSE) {
  if(is.numeric(X)) {
    X <- matrix(X, ncol=1)
  } else {
    stop("X must be a numeric vector or matrix")
  }
}

if(ncol(X)==1)
  {  X <- sort(X)
     X <- matrix(c(X[length(X)],X[-length(X)]),ncol=1)
     nouter <- 2
     outerpoints <- c(1,nrow(X))
  }
else
  {
  outerpoints <- unique(c(convhulln(X)))
  nouter <- length(outerpoints)
  }

if (is.null(initial)) y <- initialy(X)
else y <- initial
Xlong <- rbind(X[,],X[outerpoints,])
ylong <- c(y,rep(min(y)-1,nouter))
opts <- rep(0,13)
opts[1] <- as.double(-c) #-for minimization; c for initial step length
opts[2] <- ytol #distance
opts[3] <- sigmatol #function
opts[4] <- 15000 #maximum number of iterations
opts[5] <- as.double(verbose) #display control
opts[6] <- integraltol #integral
opts[7] <- as.double(alpha) #dilation factor
opts[8] <- 1*10^(-11) #for numerical gradient approx

parameters <- c(stepscale,stepscale2,stepscale3,stepscale4,desiredsize)

out <- .C("logconesteff",yvalue=as.double(y),as.double(X),as.integer(ncol(X)),as.integer(nrow(X)),options= as.double(opts),minvalue=double(1),parameters=as.double(parameters),nouter=as.integer(nouter),PACKAGE="LogConcDEAD")

y <- out$yvalue[1:nrow(X)]
Xy <- cbind(X,y)
Xymin <- cbind(X, min(y)-1)
AllPoints <- rbind(Xy, Xymin)
n <- nrow(X)
d <- ncol(X)
ConvHull <- .Call("convhulln",p=AllPoints,options="",PACKAGE="LogConcDEAD")
ConvHull <- ConvHull[apply(ConvHull<=n,1,sum)==d+1,,drop=FALSE]
nrows <- nrow(ConvHull)
beta <- NULL
b <- NULL

verts <- array(dim=c(nrows,d,d))
vertsoffset <- array(dim=c(nrows,d))
for (j in 1:nrows) {
Vertices <- Xy[ConvHull[j,],]
A <- Vertices[-1,] - rep(1,d) %o% Vertices[1,]
z <- A[,d+1]
A <- (A[,-(d+1)])
a <- Vertices[1,-(d+1)]
b <- rbind(b,solve((A),z))
beta <- rbind(beta,a%*%b[j,] - Vertices[1,d+1])
verts[j,,] <- solve(t(X[ConvHull[j,-1],])-X[ConvHull[j,1],])
vertsoffset[j,] <- verts[j,,]%*%X[ConvHull[j,1],]
}


r <- list(x=X,logMLE=out$yvalue,NumberOfEvaluations=out$options[9:11],MinSigma=out$minvalue,b=b,beta=beta,chull=ConvHull,verts=verts,vertsoffset=vertsoffset,NumberInHull=nouter)

class(r) <- "LogConcDEAD"

return(r)
}

