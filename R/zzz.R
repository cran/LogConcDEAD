

.First.lib <- function(lib, pkg) {
  if(version$major==0 && version$minor < 62)
    stop("This version for R 0.62 or later")
##  if(!require(geometry))
  ##  stop("You need the geometry package")
  library.dynam("LogConcDEAD", pkg, lib)
}
