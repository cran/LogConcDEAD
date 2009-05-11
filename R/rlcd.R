## Retain lcd.sample for consistency with old version
'lcd.sample' <- function(lcd, nsample=1) {
  warning("lcd.sample is deprecated.  Use rlcd instead")
  return( rlcd( nsample, lcd ) )
}

## For consistency with R naming convention
'rlcd' <- function(n=1, lcd) {
  triang <- lcd$triang
  x <- lcd$x
  logMLE <- lcd$logMLE
  nrows <- nrow(triang)
  d <- ncol(x)
  prob <- rep( 0, nrows )
  for (j in 1:nrows) {
    prob[j] <- lcd$detA[j]*JAD( lcd$logMLE[ triang[ j, ] ] )
  }
  prob <- prob/sum(prob)
  
  samples <- matrix( 0, nrow=n, ncol=d )
  
  ## pick a simplex for each sample
  simp <- sample( 1:nrows, n, prob = prob, replace = TRUE )
  for ( i in 1:n ) {
    while( sum( samples[ i, ] == 0 ) ) {
       ## generate point on unit simplex
      w <- rexp( d+1 )
      w <- w/sum( w )
      y <- logMLE[ triang[ simp[ i ], ] ]
      ##evaluate at the corresponding point
      fw <- exp( y %*% w )
      maxfx <- max( exp( y ) )
      u <- runif( 1 )
      if ( u < fw/maxfx ) {
        samples[ i, ] <- w%*%x[ triang[ simp[ i ], ] , ]
      }
    }
   }
  return( samples )
}














