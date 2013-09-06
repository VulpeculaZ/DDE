library(CollocInfer)
## Simulating data
vectorD.Gen <- function(t, y, parms){
    if(t<0)
        lag <- 0.05
    else
        lag <- lagvalue(t - vectorD.tau)
    dy <- vectorD.b * lag * (1 - y) - vectorD.a* y
    list(dy, dy = dy)
}


vectorD.a <- 1
vectorD.b <- 1.8
vectorD.tau <- 0.8
yinit <- 0.05
times <- seq(-0.8, 25, by = 0.1)
## solve the model

yout <- dede(y = yinit, times = times, func = vectorD.Gen, parms = NULL, atol = 1e-10)
plot(yout, which = 1, type = "l", lwd = 2, main = "Vector Disease Model")
plot(yout[,2], yout[,3], xlab = "y", ylab = "dy", type = "l", lwd = 2)

times <- times
xout <- yout[,2] + rnorm(259, sd = 0.01)
xout[xout < 0] <- 0
points(times, xout)

xout.d <- xout[times >= 0]

xout <- xout[times >=0.8]
times <- knots <- times[times >= 0.8]
norder = 3
nbasis = length(knots) + norder - 2
range  = range(knots)

vectorPars <- c(1,1.8,0.8)
names(vectorPars) <- c("a","b","tau")

vectorBasis <- create.bspline.basis(range=range,nbasis=nbasis, norder=norder,breaks=knots)

fdnames=list(NULL,c('S'),NULL)
bfdPar = fdPar(vectorBasis,lambda=1,int2Lfd(1))
DEfd0 = smooth.basis(knots,xout,bfdPar,fdnames=fdnames)$fd

par(ask=FALSE)
plotfit.fd(xout, times, DEfd0)

## extract the coefficients and assign variable names
coefs0 = DEfd0$coefs
colnames(coefs0) = "S"

##  Set a value for lambda
lambda = 1000

## Data
vectorData <- matrix(xout, ncol =1, dimnames = list(NULL,"x"))

## Setting initial values
vectorPars <- c(1,1.8,0.8)
names(vectorPars) <- c("a","b","tau")


dde.1fit <- Smooth.LS.DDE(vectorD.Fun, vectorData, times, vectorPars, coefs0, vectorBasis, lambda, in.meth='nlminb')

data <- ProfileSSE(pars, allpars, times, data, coefs, lik, proc, in.meth = "nlminb",
                   control.in = NULL, active = 1:length(pars), dcdp = NULL,
                   oldpars = NULL, use.nls = TRUE, sgn = 1)


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
