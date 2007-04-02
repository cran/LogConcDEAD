"convhulln" <-
function (p, options = " ") 
.Call("convhullnew", as.matrix(p), as.character(options),PACKAGE="LogConcDEAD")
