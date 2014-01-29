source("toload.R")
source("DSIRfn.R")
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

DSIR.pars <- c(2, 0.25, 0.5, 2, 1)
names(DSIR.pars) <- c("beta", "b","gamma", "f", "tau")
times <- seq(-DSIR.pars["tau"], 25, by = 0.1)
yinit <- c(0.85, 0.15)
lambda <- 1000
sdy <- 0.02
n <- 100

yout <- dede(y = yinit, times = times, func = DSIR.gen, parms = DSIR.pars, atol = 1e-10)

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


coefs.list <- conv.list <- data.list <- list()
par.mat <- c()
set.seed(42)
for (i in 1:n){
    initPars <- DSIR.pars[c(1, 3,5)] + runif(3, -0.2, 0.2)
    xout <- c()
    xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = sdy))
    xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = sdy))
    xout[xout < 0] <- 0
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
    dde.fit <- Profile.LS.DDE(DSIRfn, dsirData, times.d, initPars, coefs = coefs.d, basisvals = basis.d, lambda, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, tauIndex = c(FALSE, TRUE))
    DEfd2 = fd(dde.fit$coefs,basis.d, fdnames)
    coefs.list[[i]] <- dde.fit$coefs
    par.mat <- rbind(par.mat, c(dde.fit$pars, initPars))
    conv.list[[i]] <- dde.fit$convInfo
    data.list[[i]] <- dde.fit$data
}

curSeed <- .Random.seed

save(curSeed, file = "curSeed.Rdata")
save(DSIR.pars, coefs.list, conv.list, data.list, par.mat, file = paste(paste(n,sdy,lambda, sep="_"),".RData"))
