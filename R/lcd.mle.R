'lcd.mle' <- function(X, weights=rep(1/nrow(X),nrow(X)), y=initialy(X), verbose=-1, alpha=5, c=1, sigmatol=10^-8, integraltol=10^-4, ytol=10^-4, stepscale=5.1, stepscale2=2,stepscale3=1.5,stepscale4=1.05,desiredsize=3.3,Jtol=0.001) {
##  if(!require("geometry",quietly=TRUE)) stop("you need to install the geometry package")

  if (is.matrix(X)==FALSE) {
    if(is.numeric(X)) {
      X <- matrix(X, ncol=1)
    }
    else {
      stop("X must be a numeric vector or matrix")
    }
  }
  if(length(weights)!=nrow(X)) stop("there must be one w_i for each observation")
  if(sum(weights <= 0)) stop("all weights must be strictly positive")
  weights <- weights/sum(weights)
  
  if(ncol(X)==1)
    {  outerpoints <- c(which.min(X),which.max(X))
     }
  else
    {  chull <- convhullnew(X)
       outerpoints <- unique(c(chull))
     }

  innerpoints <- (1:nrow(X))[-outerpoints]
  nouter <- length(outerpoints)
  lcdsort <- c(outerpoints,innerpoints)
  
  ##Sort out the weights if necessary
  if(sum(weights <= 0)) stop("all weights must be strictly positive")
  weights <- weights/sum(weights)
  
  ##Make the initial vector
  
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

  out <- .C("logconestw",yvalue=as.double(y[lcdsort]),as.double(X[lcdsort,]),as.integer(ncol(X)),as.integer(nrow(X)),as.double(weights[lcdsort]),options= as.double(opts),minvalue=double(1),parameters=as.double(parameters),Jtol=as.double(Jtol),nouter=as.integer(nouter),PACKAGE="LogConcDEAD")

  y[lcdsort] <- out$yvalue
  Xy <- cbind(X,y)
  Xymin <- cbind(X, min(y)-1)
  AllPoints <- rbind(Xy, Xymin)
  n <- nrow(X)
  d <- ncol(X)
  ConvHull <- convhullnew(AllPoints)
  ConvHull <- ConvHull[apply(ConvHull<=n,1,sum)==d+1,,drop=FALSE]
  nrows <- nrow(ConvHull)
  beta <- NULL
  b <- NULL
  verts <- array(dim=c(nrows,d,d))
  vertsoffset <- array(dim=c(nrows,d))
  dropme <- NULL
  for (j in 1:nrows) {
    Vertices <- Xy[ConvHull[j,],]
    At <- Vertices[-1,] - rep(1,d) %o% Vertices[1,]
    z <- At[,d+1]
    At <- (At[,-(d+1),drop=FALSE])
    if(abs(det(At)) < 10^-8) {
      dropme <- c(dropme,j)
    }
    else{
    a <- Vertices[1,-(d+1)]
    bnew <- solve((At),z)
    b <- rbind(b,bnew)
    beta <- rbind(beta,a%*%bnew - Vertices[1,d+1])
    verts[j,,] <- solve(t(X[ConvHull[j,-1],])-X[ConvHull[j,1],])
    vertsoffset[j,] <- verts[j,,]%*%X[ConvHull[j,1],]
  }
  }

  rownames(b) <- NULL
  if(!is.null(dropme))  ConvHull <- ConvHull[-dropme,]


  if(ncol(X)==1) {
    midpoint = mean(X)
    chull <- c(which.min(X),which.max(X))
    outoffset <- NULL
    outnorm <- NULL
  }

  else {
    outnorm <- NULL
    outoffset <- X[chull[,1],]
    midpoint <- apply(X,2,mean) ##find a midpoint
    for (i in 1:nrow(chull)) {
      tmp <- t(X[chull[i,-1],,drop=FALSE] - rep(1,ncol(X)-1)%*%X[chull[i,1],,drop=FALSE])
      tmp <- Null(tmp)
      tmp <- -c(sign((midpoint - outoffset[i,])%*%tmp))*tmp
      outnorm <- rbind(outnorm,t(tmp))
    } 
  }
  r <- list(x=X,w=weights,logMLE=y,NumberOfEvaluations=out$options[9:11],MinSigma=out$minvalue,b=b,beta=beta,triang=ConvHull,verts=verts,vertsoffset=vertsoffset,chull=chull,outnorm=outnorm,outoffset=outoffset)
  class(r) <- "LogConcDEAD"
  return(r)}


