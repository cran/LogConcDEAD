##Note this has identical functionality to the function "convhulln" in the package "geometry"
"convhullnew" <-
function (p, options = "Qt" )
  ## Default triangulated output
.Call("convhull", as.matrix(p), as.character(options),PACKAGE="LogConcDEAD")
