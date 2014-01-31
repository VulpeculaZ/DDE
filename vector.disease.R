source("./sources/toload.R")
source("./sources/vector.D.funcs.R")
nls.control(warnOnly = TRUE)


## Function to simulating data
vectorD.Gen <- function(t, y, parms){
    if(t<0)
        lag <- 0.05
    else
        lag <- lagvalue(t - parms["tau"])
    dy <- (parms["b"] + sin(t)) * lag * (1 - y) - parms["a"]* y
    list(dy, dy = dy)
}
vectorPars <- c(1, 1.8,0.8)
names(vectorPars) <- c("a","b","tau")
yinit <- 0.05
times <- seq(-0.8, 25, by = 0.1)
## solve the model

yout <- dede(y = yinit, times = times, func = vectorD.Gen, parms = vectorPars, atol = 1e-10)
# plot(yout, which = 1, type = "l", lwd = 2, main = "Vector Disease Model")
# plot(yout[,2], yout[,3], xlab = "y", ylab = "dy", type = "l", lwd = 2)


times0 <- knots0 <- times[times >= 0]
times.d <- knots.d <- times[times >= 5]
norder = 3
nbasis.d = length(knots.d) + norder - 2
nbasis0 <- length(knots0) + norder - 2
range0  = range(knots0)
range.d <- range(knots.d)

## vectorPars <- c(1,1.8,0.8)
## names(vectorPars) <- c("a","b","tau")
vectorBasis0 <- create.bspline.basis(range=range(knots0), nbasis=nbasis0, norder=norder, breaks=knots0)
vectorBasis.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)

##  Set a value for lambda
lambda = 1000
fdnames=list(NULL,c('S'),NULL)
bfdPar0 = fdPar(vectorBasis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(vectorBasis.d,lambda=1,int2Lfd(1))

## Setting initial values
vectorPars <- c(0.8, 1.5,0.9)
names(vectorPars) <- c("a","b","tau")


xout <- yout[,2] + rnorm(length(yout[,2]), sd = 0.01)
xout[xout < 0] <- 0
xout0 <- xout[times >= 0]
xout.d <- xout[times >= 5]
DEfd0 <- smooth.basis(knots0,xout0,bfdPar0,fdnames=fdnames)$fd
DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)$fd
## temp.fit <- eval.fd(times.d, DEfd.d)
## par(ask=FALSE)
## plotfit.fd(xout[times >=0], times0, DEfd0)
## plotfit.fd(xout[times >=5], times.d, DEfd.d)
## extract the coefficients and assign variable names
coefs0 <- DEfd0$coefs
colnames(coefs0) = "S"
coefs.d <- DEfd.d$coefs
colnames(coefs.d) = "S"
## Data
vectorData <- matrix(xout.d, ncol =1, dimnames = list(NULL,"S"))


## debug(Smooth.LS.DDE)
dde.1fit2 <- Smooth.LS.DDE(vectorFun, vectorData, times.d, vectorPars, coefs = coefs.d, coefs0, basisvals = vectorBasis.d, basisvals0 =  vectorBasis0, lambda, in.meth='nlminb', delay = delay)
# Let's have a look at this
coefs1 = dde.1fit2$coefs
DEfd1 = fd(coefs1, vectorBasis.d)
plotfit.fd(vectorData,times.d , DEfd1)

## debug(Profile.LS.DDE)
# Let's have a look at this

dde.fit <- Profile.LS.DDE(vectorFun, vectorData, times.d, vectorPars, coefs = coefs1, basisvals = vectorBasis.d,  lambda, in.meth='nlminb', delay = delay, basisvals0 = vectorBasis0, coefs0 = coefs0, tauIndex = c(FALSE, FALSE, TRUE))
DEfd2 = fd(dde.fit$coefs,vectorBasis.d, fdnames)
plotfit.fd(vectorData,times.d , DEfd2)





