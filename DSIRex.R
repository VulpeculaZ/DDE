source("DSIRfn.R")

## Function to simulate  SIR model:
SIR.gen <- function(t, y, parms){
    sint <- sin(t / parms["f"]) + 1
    dyS <- -parms["beta"] * y[2] * y[1] + (1 - y[1]) * (parms["b"]* sint + 0.1)
    dyI <- parms["beta"] * y[2] * y[1] - (parms["gamma"] + (parms["b"] * sint + 0.1)) * y[2]
    list(c(dyS, dyI))
}


## Function to simulate DSIR
DSIR.gen <- function(t, y, parms){
    if(t < 0){
        lagI <- 0.05
    }
    else
        lagI <- lagvalue(t - parms["tau"], 2)
    sint <- sin(t / parms["f"]) + 1
    if(y[2] > 1)
        y[2] <- 1
    if(y[2] < 0)
        y[2] <- 0
    dyS <- -parms["beta"] * lagI * y[1] + (1 - y[1]) * (parms["b"] * sint + 0.1)
    dyI <- parms["beta"] * lagI * y[1] - (parms["gamma"] + (parms["b"] * sint + 0.1)) * y[2]
    list(c(dyS, dyI))
}


yinit <- c(0.95, 0.05)
times <- seq(0, 25, by = 0.1)
## Simulation for SIR
SIR.pars <- c(2, 0.25, 0.5, 1)
names(SIR.pars) <- c("beta", "b","gamma", "f")
yout <- ode(y = yinit, times = times, func = SIR.gen, parms = SIR.pars, atol = 1e-10)
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, main = "SIR Model")

## Simulation for DSIR
tau <- 1
DSIR.pars <- c(SIR.pars, tau)
names(DSIR.pars) <- c("beta", "b","gamma", "f", "tau")
times <- seq(-DSIR.pars["tau"], max(times), by = 0.1)
yout <- dede(y = yinit, times = times, func = DSIR.gen, parms = DSIR.pars, atol = 1e-10)
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, main = "Delayed SIR Model")


times0 <- knots0 <- times[times >= 0]
times.d <- knots.d <- times[times >= 5]
norder = 3
nbasis.d = length(knots.d) + norder - 2
nbasis0 <- length(knots0) + norder - 2
range0  = range(knots0)
range.d <- range(knots.d)
basis0 <- create.bspline.basis(range=range(knots0), nbasis=nbasis0, norder=norder, breaks=knots0)
basis.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar0 = fdPar(basis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(basis.d,lambda=1,int2Lfd(1))


xout <- c()
xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 0.05))
xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 0.05))
## points(times, xout)
xout0 <- xout[times >= 0,]
xout.d <- xout[times >= 5,]
DEfd0 <- smooth.basis(knots0,xout0,bfdPar0,fdnames=fdnames)$fd
DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)$fd
## temp.fit <- eval.fd(times.d, DEfd.d)
## par(ask=FALSE)
plotfit.fd(xout[times >=0,], times0, DEfd0)
plotfit.fd(xout[times >=5,], times.d, DEfd.d)

## extract the coefficients and assign variable names
coefs0 <- DEfd0$coefs
colnames(coefs0) = c("S","I")
coefs.d <- DEfd.d$coefs
colnames(coefs.d) = c("S", "I")


##  Set a value for lambda
lambda = 1000

## Data
dsirData <- matrix(xout.d, ncol =2, dimnames = list(NULL,c("S", "I")))

## Setting initial values
initPars <- c(1.4, 0.6, 0.9)
names(initPars) <- c("beta","gamma","tau")

# nls.control(warnOnly = TRUE)
# debug(Smooth.LS.DDE)
dde.1fit2 <- Smooth.LS.DDE(DSIRfn, dsirData, times.d, initPars, coefs = coefs.d, coefs0, basisvals = basis.d, basisvals0 =  basis0, lambda, in.meth='nlminb',tauMax = 5, delay = delay)
# Let's have a look at this

coefs1 = dde.1fit2$coefs
DEfd1 = fd(coefs1, basis.d)
plotfit.fd(dsirData, times.d , DEfd1)

dde.fit <- Profile.LS.DDE(DSIRfn, dsirData, times.d, initPars, coefs = coefs1, basisvals = basis.d, lambda, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, tauIndex = c(FALSE, TRUE))

DEfd2 = fd(dde.fit$coefs,basis.d, fdnames)
plotfit.fd(dsirData, times.d , DEfd2)

fit.list <- list()
par.mat <- c()
set.seed(42)
for (i in 1:20){
    xout <- c()
    xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 0.01))
    xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 0.01))
    xout0 <- xout[times >= 0,]
    xout.d <- xout[times >= 5,]
    DEfd0 <- smooth.basis(knots0,xout0,bfdPar0,fdnames=fdnames)$fd
    DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)$fd
    ## extract the coefficients and assign variable names
    coefs0 <- DEfd0$coefs
    colnames(coefs0) = c("S","I")
    coefs.d <- DEfd.d$coefs
    colnames(coefs.d) = c("S", "I")
    dsirData <- matrix(xout.d, ncol =2, dimnames = list(NULL,c("S", "I")))
    dde.fit <- Profile.LS.DDE(DSIRfn, dsirData, times.d, initPars, coefs = coefs0, basisvals = basis.d, lambda, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, tauIndex = c(FALSE, TRUE))
    DEfd2 = fd(dde.fit$coefs,basis.d, fdnames)
    fit.list[[i]] <- DEfd2
    par.mat <- rbind(par.mat, dde.fit$pars)
}

hist(par.mat[,1])
hist(par.mat[,2])
hist(par.mat[,3])
## Varying initial values, 1000 replicas
