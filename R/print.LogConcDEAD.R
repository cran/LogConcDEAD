print.LogConcDEAD <- function(x, ...)
  {
niter <- x$NumberOfEvaluations[1]
if (niter==-1) errormessage <- "allocation error"
else if (niter==-2) errormessage <- "improper space dimension"
else if (niter==-3) errormessage <- "sigma returns an improper value"
else if (niter==-4) errormessage <- "grad_sigma returns a zero or improper vector at the starting point"
else if(niter==-7) errormessage <- "sigma is unbounded"
else if (niter==-8) errormessage <- "gradient is zero, but stopping criteria are not fulfilled"
else if (niter==-9) errormessage <- "iterations limit exceeded"
else if (niter==-11) errormessage <- "Premature stop is possible"
else if (niter==-12) errormessage <- "Result may not provide the true optimum"
else if (niter==-13) errormessage <- "Result may be inaccurate in view of a point"
else if (niter==-14) errormessage <- "Result may be inaccurate in view of a function value"

if(niter>0) cat("\n Log MLE at observations: \n", (x$logMLE),
    "\n\n Number of Iterations: ",x$NumberOfEvaluations[1],
    "\n\n Number of Function Evaluations: ",x$NumberOfEvaluations[2],
    fill=10)
else cat("SolvOpt error: ",errormessage, "\n")
}
