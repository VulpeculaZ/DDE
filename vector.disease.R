library(CollocInfer)
library(deSolve)

## Function to simulating data
vectorD.Gen <- function(t, y, parms){
    if(t<0)
        lag <- 0.05
    else
        lag <- lagvalue(t - vectorD.tau)
    dy <- (vectorD.b + sin(t)) * lag * (1 - y) - vectorD.a* y
    list(dy, dy = dy)
}
vectorD.a <- 1
vectorD.b <- 1.8
vectorD.tau <- 0.8
yinit <- 0.05
times <- seq(-0.8, 50, by = 0.1)
## solve the model

yout <- dede(y = yinit, times = times, func = vectorD.Gen, parms = NULL, atol = 1e-10)
plot(yout, which = 1, type = "l", lwd = 2, main = "Vector Disease Model")
plot(yout[,2], yout[,3], xlab = "y", ylab = "dy", type = "l", lwd = 2)

xout <- yout[,2] + rnorm(length(yout[,2]), sd = 0.005)
xout[xout < 0] <- 0
points(times, xout)
xout0 <- xout[times >= 0]

xout.d <- xout[times >= 5]
times0 <- knots0 <- times[times >= 0]
times.d <- knots.d <- times[times >= 5]
norder = 3
nbasis.d = length(knots.d) + norder - 2
nbasis0 <- length(knots0) + norder - 2
range0  = range(knots0)
range.d <- range(knots.d)

vectorPars <- c(1,1.8,0.8)
names(vectorPars) <- c("a","b","tau")
vectorBasis0 <- create.bspline.basis(range=range(knots0), nbasis=nbasis0, norder=norder, breaks=knots0)
vectorBasis.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)


fdnames=list(NULL,c('S'),NULL)
bfdPar0 = fdPar(vectorBasis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(vectorBasis.d,lambda=1,int2Lfd(1))
DEfd0 = smooth.basis(knots0,xout0,bfdPar0,fdnames=fdnames)$fd
DEfd0.true <- smooth.basis(knots0,yout[times >=0,2],bfdPar0,fdnames=fdnames)$fd
DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)$fd
temp.fit <- eval.fd(times.d, DEfd.d)
par(ask=FALSE)
plotfit.fd(yout[times >=0,2], times0, DEfd0.true)
plotfit.fd(xout[times >=5], times.d, DEfd.d)

## extract the coefficients and assign variable names
coefs0 <- DEfd0.true$coefs
colnames(coefs0) = "S"
coefs.d <- DEfd.d$coefs
colnames(coefs.d) = "S"

##  Set a value for lambda
lambda = 1000

## Data
vectorData <- matrix(xout.d, ncol =1, dimnames = list(NULL,"S"))

## Setting initial values
vectorPars <- c(0.8, 1.5,0.8)
names(vectorPars) <- c("a","b","tau")

## debug(Smooth.LS.DDE)
dde.1fit2 <- Smooth.LS.DDE(vectorFun, vectorData, times.d, vectorPars, coefs = coefs.d, coefs0, basisvals = vectorBasis.d, basisvals0 =  vectorBasis0, lambda, in.meth='nlminb', delay = delay)
# Let's have a look at this
coefs1 = dde.1fit2$coefs
DEfd1 = fd(coefs1, vectorBasis.d)
plotfit.fd(vectorData,times.d , DEfd1)

## debug(Profile.LS.DDE)
# Let's have a look at this
dde.fit <- Profile.LS.DDE(vectorFun, vectorData, times.d, vectorPars, coefs = coefs1, basisvals = vectorBasis.d,  lambda, in.meth='nlminb', delay = delay, basisvals0 = vectorBasis0, coefs0 = coefs0, tauIndex = c(FALSE, FALSE, TRUE))
print("")
print(dde.fit$pars)
DEfd2 = fd(dde.fit$coefs,vectorBasis.d, fdnames)
plotfit.fd(vectorData,times.d , DEfd2)


save(list = ls(), file = "nls1.Rdata")
## Testing a intermediate ProfileSSE.Allpar.DDE
load("nls1.Rdata")
tmpfd <- fd(ncoefs, basisvals)
plotfit.fd(data, times, tmpfd)
    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names

    ## Squared errors: No need to change for DDE
    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
    f = as.vector( as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more))*sqrt(weights))
    isnaf = is.na(f)
    f[isnaf] = 0

    dlikdp = lik$more$dfdp(times,devals,pars,lik$more$more)
    dlikdp = matrix(dlikdp,dim(dlikdp)[1]*dim(dlikdp)[2],dim(dlikdp)[3])
    ## dlikdp will be zero if the likelihood doesn't directly have ode parameters,
    ## which is true for least square case.
    dlikdx = lik$more$dfdx(times,devals,pars,lik$more$more)
    ## dlikdx[i,,] is an identity matrix for every i.

    dlikdc = c()
    for(i in 1:dim(dlikdx)[2]){
        tH = c()
        for(j in 1:dim(dlikdx)[3]){
            ## Only dlikdx[,i,i] are non-zero (all 1's)
            tH = cbind(tH,as.matrix(diag(dlikdx[,i,j])%*%lik$bvals))
        }
        dlikdc = rbind(dlikdc,tH)
    }

    ## ??dlikdc: why 0.5 and 0.5 ??
    d2Hdc2  = SplineCoefsDC2sparse(ncoefs,times,data,lik,proc,pars)
    ## need not be changed?
    d2Hdcdp = SplineCoefsDCDP(ncoefs,times,data,lik,proc,pars)
    ## Got warning message:
    ## In dim(weights[whichrows, ]) == dim(diffs):
    ## longer object length is not a multiple of shorter object length?

    ## Use Implicit function theorem:
    if(is.matrix(d2Hdc2)){
        ## When will it not be a matrix? How to use solve in that case?
        dcdp = ginv(d2Hdc2) %*% d2Hdcdp
    } else {
        dcdp = as.matrix(solve(d2Hdc2,d2Hdcdp))
    }

    ## Chain rule:
    df = dlikdc%*%dcdp + dlikdp
    df[isnaf,] = 0
    colnames(df) = proc$more$parnames





daybasis <- create.fourier.basis(c(0, 365), nbasis=65)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#  See analyses of weather data.
#  Set up sampling points at mid days
#  Convert the data to a functional data object
tempfd <- smooth.basis(day.5,  CanadianWeather$dailyAv[,,"Temperature.C"],
                       daybasis)$fd
#   set up the harmonic acceleration operator
Lbasis  <- create.constant.basis(c(0, 365))
Lcoef   <- matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj  <- fd(Lcoef,Lbasis)
bwtlist <- fd2list(bfdobj)
harmaccelLfd <- Lfd(3, bwtlist)
#   evaluate the value of the harmonic acceleration
#   operator at the sampling points
Ltempmat <- eval.fd(day.5, tempfd, harmaccelLfd)



proc$d2fdc2 <- function(coefs, proc$bvals, pars, proc$more)
{
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    H1 = make.SSElik()$d2fdx2(ddevals, more$qpts, devals, pars,
        more)
    H2 = more$dfdx(more$qpts, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, H1[,
        , 1, drop = TRUE])
    H = list(len = dim(bvals$bvals)[2])
    for (i in 1:dim(devals)[2]) {
        H[[i]] = list(len = dim(devals))
        for (j in 1:dim(devals)[2]) {
            H[[i]][[j]] = t(bvals$bvals) %*% diag(H1[, i, j]) %*%
                bvals$bvals - 2 * t(bvals$dbvals) %*% diag(H2[,
                i, j] * weights[, i]) %*% bvals$bvals - 2 * t(bvals$bvals) %*%
                diag(H2[, j, i] * weights[, j]) %*% bvals$dbvals
        }
        H[[i]][[i]] = H[[i]][[i]] + 2 * t(bvals$dbvals) %*% diag(weights[,
            i]) %*% bvals$dbvals
    }
    H = blocks2mat(H)
    return(H)
}

temp.d2dfdx2 <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx = more$dfdx(times, devals, pars, more$more)
    d2fdx2 = more$d2fdx2(times, devals, pars, more$more)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdx2)[3]) {
        for (j in 1:dim(d2fdx2)[4]) {
            H[, i, j] = apply(-difs * d2fdx2[, , i, j] + weights *
             dfdx[, , j] * dfdx[, , i], 1, sum)
        }
    }
    return(2 * H)
}

## debugging in: SplineCoefsDC2sparse(coefs, times, data, lik, proc, pars, sgn):
debug(proc$d2fdc2)
temp <- proc$d2fdc2(coefs2, proc$bvals, pars, proc$more)

##
debug(temp.d2dfdx2)
H1 = temp.d2dfdx2(ddevals, more$qpts, devals, pars, more)

collocNS <- getNamespace("CollocInfer")
unlockBinding("checkweights", collocNS)
assignInNamespace("checkweights",checkweights_1d, ns="CollocInfer", envir=collocNS)
assign("checkweights", checkweights_1d, envir=collocNS)
lockBinding("checkweights",collocNS)
ckeckweights <- checkweights_1d
rm() ## get rid of global copy
