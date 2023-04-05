### R code from vignette source 'LogConcDEAD.Rnw'

###################################################
### code chunk number 1: LogConcDEAD.Rnw:126-130
###################################################
options(prompt="R> ")
require( "LogConcDEAD" )
require( "mvtnorm" )
options(width=72)


###################################################
### code chunk number 2: setn
###################################################
n <- 100


###################################################
### code chunk number 3: LogConcDEAD.Rnw:416-420
###################################################
set.seed(1)
n <- 100
x <- sort(rgamma(n,shape=2))
out1 <- mlelcd(x) ## LogConcDEAD estimate


###################################################
### code chunk number 4: plot:1d (eval = FALSE)
###################################################
## ylim <- c(0,0.4) 
## lgdtxt <- c("LogConcDEAD", "true")
## lgdlty <- c(1,3)
## plot(out1, ylim=ylim,lty=1)
## lines(x, x*exp(-x), lty=3)
## legend(x=3, y=0.4, lgdtxt, lty=lgdlty)


###################################################
### code chunk number 5: plot:log1d (eval = FALSE)
###################################################
## ylim <- c(-4,-1)
## lgdtxt <- c("LogConcDEAD",  "true")
## lgdlty <- c(1,3)
## plot(out1, uselog=TRUE, lty=1)
## lines(x, log(x)-x, lty=3)
## legend(x=3, y=-1, lgdtxt, lty=lgdlty) 


###################################################
### code chunk number 6: fig_1d
###################################################
ylim <- c(0,0.4) 
lgdtxt <- c("LogConcDEAD", "true")
lgdlty <- c(1,3)
plot(out1, ylim=ylim,lty=1)
lines(x, x*exp(-x), lty=3)
legend(x=3, y=0.4, lgdtxt, lty=lgdlty)


###################################################
### code chunk number 7: fig_log1d
###################################################
ylim <- c(-4,-1)
lgdtxt <- c("LogConcDEAD",  "true")
lgdlty <- c(1,3)
plot(out1, uselog=TRUE, lty=1)
lines(x, log(x)-x, lty=3)
legend(x=3, y=-1, lgdtxt, lty=lgdlty) 


###################################################
### code chunk number 8: setn2
###################################################
n <- 100


###################################################
### code chunk number 9: LogConcDEAD.Rnw:487-491
###################################################
set.seed(22)
d <- 2
n <- 100
x <- matrix(rnorm(n*d),ncol=d)


###################################################
### code chunk number 10: LogConcDEAD.Rnw:503-504
###################################################
out <- mlelcd(x,verbose=50)


###################################################
### code chunk number 11: LogConcDEAD.Rnw:550-551
###################################################
names(out)


###################################################
### code chunk number 12: LogConcDEAD.Rnw:559-560
###################################################
out$logMLE[1:5]


###################################################
### code chunk number 13: LogConcDEAD.Rnw:572-573
###################################################
out$triang[1:5,]


###################################################
### code chunk number 14: LogConcDEAD.Rnw:580-582
###################################################
out$b[1:5,]
out$beta[1:5]


###################################################
### code chunk number 15: LogConcDEAD.Rnw:601-603
###################################################
out$NumberOfEvaluations
out$MinSigma


###################################################
### code chunk number 16: LogConcDEAD.Rnw:607-608
###################################################
out$chull[1:5,]


###################################################
### code chunk number 17: LogConcDEAD.Rnw:614-616
###################################################
out$outnorm[1:5,]
out$outoffset[1:5,]


###################################################
### code chunk number 18: LogConcDEAD.Rnw:646-649
###################################################
g <- interplcd(out, gridlen=200)
g1 <- interpmarglcd(out, marg=1)
g2 <- interpmarglcd(out, marg=2)


###################################################
### code chunk number 19: plot:2d (eval = FALSE)
###################################################
## par(mfrow=c(1,2), pty="s", cex=0.7) #square plots
## plot(out,g=g,addp=FALSE,asp=1)
## plot(out,g=g,uselog=TRUE,addp=FALSE,asp=1)


###################################################
### code chunk number 20: LogConcDEAD.Rnw:666-671
###################################################
png(file="LogConcDEAD-fig_2d.png")
oldpar <- par(mfrow = c(1,1))
par(mfrow=c(1,2), pty="s", cex=0.7) #square plots
plot(out,g=g,addp=FALSE,asp=1)
plot(out,g=g,uselog=TRUE,addp=FALSE,asp=1)
par(oldpar)
dev.off()


###################################################
### code chunk number 21: plot:2dmarg (eval = FALSE)
###################################################
## par(mfrow=c(1,2), pty="s", cex=0.7) #normal proportions 
## plot(out,marg=1,g.marg=g1)
## plot(out,marg=2,g.marg=g2)


###################################################
### code chunk number 22: fig_2dmarg
###################################################
oldpar <- par(mfrow = c(1,1))
par(mfrow=c(1,2), pty="s", cex=0.7) #normal proportions 
plot(out,marg=1,g.marg=g1)
plot(out,marg=2,g.marg=g2)
par(oldpar)


###################################################
### code chunk number 23: plot:rgl (eval = FALSE)
###################################################
## plot(out,g=g,type="r")


###################################################
### code chunk number 24: plot:rgllog (eval = FALSE)
###################################################
## plot(out,g=g,type="r",uselog=TRUE)


###################################################
### code chunk number 25: LogConcDEAD.Rnw:732-735 (eval = FALSE)
###################################################
## plot(out,g=g,type="r")
## par3d(windowRect = c(55,66,311+256, 322+256))
## rgl.snapshot(file="rglfig.png")


###################################################
### code chunk number 26: LogConcDEAD.Rnw:744-747 (eval = FALSE)
###################################################
## plot(out,g=g,type="r",uselog=TRUE)
## par3d(windowRect = c(55,66,311+256, 322+256))
## rgl.snapshot(file="rgllog.png")


###################################################
### code chunk number 27: LogConcDEAD.Rnw:774-776
###################################################
nsamp <- 1000
mysamp <- rlcd(nsamp,out)


###################################################
### code chunk number 28: LogConcDEAD.Rnw:782-784
###################################################
colMeans(mysamp)
cov(mysamp)


###################################################
### code chunk number 29: LogConcDEAD.Rnw:791-793
###################################################
m <- 10
mypoints <- 1.5*matrix(rnorm(m*d),ncol=d)


###################################################
### code chunk number 30: LogConcDEAD.Rnw:797-798
###################################################
dlcd(mypoints,out)


###################################################
### code chunk number 31: LogConcDEAD.Rnw:813-817
###################################################
myval <- sort(dlcd(mysamp,out))
alpha <- c(.25,.5,.75)
myval[(1-alpha)*nsamp]



###################################################
### code chunk number 32: setn3
###################################################
n <- 100


###################################################
### code chunk number 33: LogConcDEAD.Rnw:832-834 (eval = FALSE)
###################################################
## install.packages("mvtnorm")
## library("mvtnorm")


###################################################
### code chunk number 34: LogConcDEAD.Rnw:837-843
###################################################
set.seed(333)
sigma <- matrix(c(1,0.2,0.2,1),nrow=2)
d <- 2
n <- 100
y <- rmvnorm(n,sigma=0.1*sigma)
xall <- round(y,digits=1)


###################################################
### code chunk number 35: LogConcDEAD.Rnw:856-859
###################################################
tmpw <- getweights(xall)
outw <- mlelcd(tmpw$x,w=tmpw$w)
gw <- interplcd(outw, gridlen=200)


###################################################
### code chunk number 36: plot:bin (eval = FALSE)
###################################################
## par(mfrow=c(1,2), pty="s", cex=0.7) #2 square plots 
## plot(outw,g=gw,asp=1,drawlabels=FALSE)
## plot(outw,g=gw,uselog=TRUE,asp=1,drawlabels=FALSE)


###################################################
### code chunk number 37: LogConcDEAD.Rnw:877-882
###################################################
png(file="LogConcDEAD-fig_bin.png")
oldpar <- par(mfrow = c(1,1))
par(mfrow=c(1,2), pty="s", cex=0.7) #2 square plots 
plot(outw,g=gw,asp=1,drawlabels=FALSE)
plot(outw,g=gw,uselog=TRUE,asp=1,drawlabels=FALSE)
par(oldpar)
dev.off()


###################################################
### code chunk number 38: setn4
###################################################
n <- 100


###################################################
### code chunk number 39: LogConcDEAD.Rnw:903-908
###################################################
set.seed(4444)
d <- 3
n <- 100
x <- matrix(rgamma(n*d,shape=2),ncol=d)
out3 <- mlelcd(x)


###################################################
### code chunk number 40: LogConcDEAD.Rnw:915-917
###################################################
mypoints <- c(0,2,4)
dmarglcd(mypoints, out3, marg=1)


###################################################
### code chunk number 41: plot:3deg (eval = FALSE)
###################################################
## par(mfrow=c(2,2),cex=0.8)
## plot(out3, marg=1)
## plot(out3, marg=2)
## plot(out3, marg=3)
## tmp <- seq(min(out3$x), max(out3$x),len=100)
## plot(tmp, dgamma(tmp,shape=2), type="l", 
## xlab="X", ylab="true marginal density")
## title(main="True density")


###################################################
### code chunk number 42: fig_3deg
###################################################
oldpar <- par(mfrow = c(1,1))
par(mfrow=c(2,2),cex=0.8)
plot(out3, marg=1)
plot(out3, marg=2)
plot(out3, marg=3)
tmp <- seq(min(out3$x), max(out3$x),len=100)
plot(tmp, dgamma(tmp,shape=2), type="l", 
xlab="X", ylab="true marginal density")
title(main="True density")
par(oldpar)


