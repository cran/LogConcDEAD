##Note this has identical functionality to the function "convhulln" in the package "geometry"
"convhullnew" <-
function (p, options = " ") 
.Call("convhullnew", as.matrix(p), as.character(options),PACKAGE="LogConcDEAD")
