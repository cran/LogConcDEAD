'initialf' <- 
function(x,n,Xin,h){
mean(apply(dnorm(rep(1,n)%o%x-Xin,sd = rep(1,n)%o%h),1,prod))
}

'initialy' <-
function(Xin){
n <- nrow(Xin)
d <- ncol(Xin)
h <- apply(Xin,2,sd)*n^(-1/(d+4))
log(apply(Xin,1,initialf,n=n,Xin=Xin,h=h))}
