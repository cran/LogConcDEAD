'mlelcd' <- function(x, w=rep(1/nrow(x),nrow(x)), y=initialy(x), verbose=-1, alpha=5, c=1, sigmatol=10^-8, integraltol=10^-4, ytol=10^-4, stepscale=5.1, stepscale2=2,stepscale3=1.5,stepscale4=1.05,desiredsize=3.3,Jtol=0.001) {
##  if(!require("geometry",quietly=TRUE)) stop("you need to install the geometry package")

  if (is.matrix(x)==FALSE) {
    if(is.numeric(x)) {
      x <- matrix(x, ncol=1)
    }
    else {
      stop("x must be a numeric vector or matrix")
    }
  }
  if(length(w)!=nrow(x)) stop("there must be one w_i for each observation")
  if(sum(w <= 0)) stop("all weights must be strictly positive")
  w <- w/sum(w)
  
  if(ncol(x)==1)
    {  outerpoints <- c(which.min(x),which.max(x))
     }
  else
    {  chull <- convhullnew(x)
       outerpoints <- unique(c(chull))
     }

  innerpoints <- (1:nrow(x))[-outerpoints]
  nouter <- length(outerpoints)
  lcdsort <- c(outerpoints,innerpoints)
  
  ##Sort out the weights if necessary
  if(sum(w <= 0)) stop("all weights must be strictly positive")
  w <- w/sum(w)
  
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

  out <- .C("logconestw",yvalue=as.double(y[lcdsort]),as.double(x[lcdsort,]),as.integer(ncol(x)),as.integer(nrow(x)),as.double(w[lcdsort]),options= as.double(opts),minvalue=double(1),parameters=as.double(parameters),Jtol=as.double(Jtol),nouter=as.integer(nouter),PACKAGE="LogConcDEAD")

  y[lcdsort] <- out$yvalue
  xy <- cbind(x,y)
  xymin <- cbind(x, min(y)-1)
  AllPoints <- rbind(xy, xymin)
  n <- nrow(x)
  d <- ncol(x)
  ConvHull <- convhullnew(AllPoints)
  ConvHull <- ConvHull[apply(ConvHull<=n,1,sum)==d+1,,drop=FALSE]
  nrows <- nrow(ConvHull)
  beta <- NULL
  b <- NULL
  verts <- array(dim=c(nrows,d,d))
  vertsoffset <- array(dim=c(nrows,d))
  dropme <- NULL
  for (j in 1:nrows) {
    Vertices <- xy[ConvHull[j,],]
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
    verts[j,,] <- solve(t(x[ConvHull[j,-1],])-x[ConvHull[j,1],])
    vertsoffset[j,] <- verts[j,,]%*%x[ConvHull[j,1],]
  }
  }

  rownames(b) <- NULL
  if(!is.null(dropme))  ConvHull <- ConvHull[-dropme,]


  if(ncol(x)==1) {
    midpoint = mean(x)
    chull <- c(which.min(x),which.max(x))
    outoffset <- matrix(c(min(x),max(x)),ncol=1)
    outnorm <- matrix(c(-1,1),ncol=1)
  }
  else {
    outnorm <- NULL
    outoffset <- x[chull[,1],]
    midpoint <- apply(x,2,mean) ##find a midpoint
    for (i in 1:nrow(chull)) {
      tmp <- t(x[chull[i,-1],,drop=FALSE] - rep(1,ncol(x)-1)%*%x[chull[i,1],,drop=FALSE])
      tmp <- Null(tmp)
      tmp <- -c(sign((midpoint - outoffset[i,])%*%tmp))*tmp
      outnorm <- rbind(outnorm,t(tmp))
    } 
  }
  r <- list(x=x,w=w,logMLE=y,NumberOfEvaluations=out$options[9:11],MinSigma=out$minvalue,b=b,beta=beta,triang=ConvHull,verts=verts,vertsoffset=vertsoffset,chull=chull,outnorm=outnorm,outoffset=outoffset)
  class(r) <- "LogConcDEAD"
  return(r)}

## retained for compatibility with previous versions
'lcd.mle' <- function(x, w=rep(1/nrow(x),nrow(x)),
y=initialy(x), verbose=-1, alpha=5, c=1, sigmatol=10^-8,
integraltol=10^-4, ytol=10^-4, stepscale=5.1,
stepscale2=2,stepscale3=1.5,stepscale4=1.05,desiredsize=3.3,Jtol=0.001)
{
warning("lcd.mle is deprecated.  Use mlelcd instead")
  return( mlelcd( x, w, y, verbose, alpha, c, sigmatol, integraltol,
ytol, stepscale, stepscale2, stepscale3 ,stepscale4 ,desiredsize, Jtol
) )
}
