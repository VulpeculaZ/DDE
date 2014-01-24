SplineCoefsDC <- function(){
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals,
        pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars,
        proc$more)
    g = as.vector(g)
    return(sgn * g)
}

proc$d2fdc2 <- function (coefs, bvals, pars, more)
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

proc$d2fdcdp <- function (coefs, bvals, pars, more)
{
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    H1 = make.SSElik()$d2fdxdp(ddevals, more$qpts, devals, pars,
        more)
    H2 = 2 * more$dfdp(more$qpts, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, H1[,
        , 1, drop = FALSE])
    H = c()
    for (i in 1:length(pars)) {
        H = cbind(H, as.vector(as.matrix(t(bvals$bvals) %*% H1[,
            , i] - t(bvals$dbvals) %*% (weights * H2[, , i]))))
    }
    return(H)
}

make.SSElik()$d2fdxdp <- function (data, times, devals, pars, more)
{
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdp = more$dfdp(times, devals, pars, more$more)
    d2fdxdp = more$d2fdxdp(times, devals, pars, more$more)
    H = array(0, c(dim(devals), length(pars)))
    for (i in 1:dim(d2fdxdp)[3]) {
        for (j in 1:dim(d2fdxdp)[4]) {
            H[, i, j] = apply(-difs * d2fdxdp[, , i, j] + weights *
                dfdx[, , i] * dfdp[, , j], 1, sum)
        }
    }
    return(2 * H)
}

SplineCoefsDCDP <- function(coefs, times, data, lik, proc, pars, sgn = 1)
{
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    d2lik = lik$d2fdxdp(data, times, devals, pars, lik$more)
    H = c()
    for (i in 1:length(pars)) {
        H = cbind(H, as.vector(as.matrix(t(lik$bvals) %*% d2lik[,
            , i])))
    }
    ## H has all 0's up to this point.
    H = H + proc$d2fdcdp(coefs2, proc$bvals, pars, proc$more)
    return(as.matrix(sgn * H))
}


lik$d2fdxdp <- function (data, times, devals, pars, more)
{
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdp = more$dfdp(times, devals, pars, more$more)
    d2fdxdp = more$d2fdxdp(times, devals, pars, more$more)
    H = array(0, c(dim(devals), length(pars)))
    for (i in 1:dim(d2fdxdp)[3]) {
        for (j in 1:dim(d2fdxdp)[4]) {
            H[, i, j] = apply(-difs * d2fdxdp[, , i, j] + weights *
                dfdx[, , i] * dfdp[, , j], 1, sum)
        }
    }
    return(2 * H)
}

 make.SSElik()$d2fdx2 <- function (data, times, devals, pars, more)
{
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

## make.SSElik()$d2fdx2 is the same as lik$d2fdx2 (in the case of LS?!)
## H[,i,i] will be 2 * identity matrix
lik$d2fdx2.DDE <- function (data, times, devals, pars, more)
{
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx = more$dfdx(times, devals, pars, more$more)
    ## d2fdx2 is a zero array
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



## Zero's
lik$d2fdxdp.d <- function (data, times, devals, pars, more)
{
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdp = more$dfdp(times, devals, pars, more$more) #0's
    d2fdxdp = more$d2fdxdp(times, devals, pars, more$more) #0's
    H = array(0, c(dim(devals), length(pars)))
    for (i in 1:dim(d2fdxdp)[3]) {
        for (j in 1:dim(d2fdxdp)[4]) {
            H[, i, j] = apply(-difs * d2fdxdp[, , i, j] + weights *
                dfdx[, , i] * dfdp[, , j], 1, sum)
        }
    }
    ## All 0's
    return(2 * H)
}

##


delay$d2fdx2(ddevals, more$qpts, devals, pars, more, devals.d, ddevals.d)



weights = checkweights(more$weights, more$whichobs, H1[,, 1, drop = TRUE])





delay$d2fdx.ddp <- function(ddevals, more$qpts, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx.d = more$dfdx.d(times, devals, pars, more$more)
    dfdp = more$dfdp(times, devals, pars, more$more)
    d2fdx.ddp = more$d2fdxdp(times, devals, pars, more$more)
    H = array(0, c(dim(devals), length(pars)))
    for (i in 1:dim(d2fdxdp)[3]) {
        for (j in 1:dim(d2fdxdp)[4]) {
            H[, i, j] = apply(-difs * d2fdx.ddp[, , i, j] + weights *
             dfdx.d[, , i] * dfdp[, , j], 1, sum)
        }
    }
    return(2 * H)
}

## Same???
delay$d2fdxdd <- function(ddevals, more$qpts, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdd = more$dfdd(times, devals, pars, more$more)
    d2fdxdd = more$d2fdxdp(times, devals, pars, more$more)
    H = array(0, c(dim(devals), length(pars)))
    for (i in 1:dim(d2fdxdp)[3]) {
        for (j in 1:dim(d2fdxdp)[4]) {
            H[, i, j] = apply(-difs * d2fdx.ddp[, , i, j] + weights *
                dfdx.d[, , i] * dfdp[, , j], 1, sum)
        }
    }
    return(2 * H)
}

delay$d2fdx.ddd <- function(ddevals, more$qpts, devals, pars, more, devals.d){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdd = more$dfdd(times, devals, pars, more$more)
    d2fdxdd = more$d2fdxdp(times, devals, pars, more$more)
    H = array(0, c(dim(devals), length(pars)))
    for (i in 1:dim(d2fdxdp)[3]){
        for (j in 1:dim(d2fdxdp)[4]) {
            H[, i, j] = apply(-difs * d2fdx.ddp[, , i, j] + weights *
                dfdx.d[, , i] * dfdp[, , j], 1, sum)
        }
    }
    return(2 * H)
}

delay$dfdd <- function(more$qpts, devals, pars, more$more, devals.d)
{
    devals = as.matrix(more$more$bvals %*% coefs)
    ddevals = as.matrix(more$more$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    g = delay$dfddMat(ddevals, more$qpts, devals, pars,
        more)
    g = apply(g, 2, sum)
    return(g)
}

delay$dfddMat <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdp = more$dfdp(times, devals, pars, more$more)
    g = c()
    for (i in 1:dim(dfdp)[3]) {
        g = cbind(g, apply(difs * dfdd[, , i], 1, sum))
    }
    return(-2 * g)
}

df.auto <- matrix(0, dim(df)[1], dim(df)[2])
eps <- 1E-4
for(i in 1:length(pars)){
    pars.i1 <- pars.i2 <- pars
    pars.i1[i] <- pars[i] + eps
    pars.i2[i] <- pars[i] - eps
    Ires.i1 <- inneropt.DDE(data,times,pars.i1,coefs,lik,proc, in.meth,control.in, basisvals = basisvals, fdobj0 = fdobj0)
    Ires.i2 <- inneropt.DDE(data,times,pars.i2,coefs,lik,proc, in.meth,control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ncoefs.i1 <- Ires.i1$coefs
    ncoefs.i2 <- Ires.i2$coefs
    devals.i1 = as.matrix(lik$bvals%*%ncoefs.i1)
    devals.i2 = as.matrix(lik$bvals%*%ncoefs.i2)
    colnames(devals.i1) <- colnames(devals.i2) <- proc$more$names
    f.i1 = as.vector(lik$more$fn(times, devals.i1, pars, lik$more$more))*sqrt(weights)
    f.i2 = as.vector(lik$more$fn(times, devals.i2, pars, lik$more$more))*sqrt(weights)
    df.auto[,i] <- (f.i2 - f.i1) / 2 / eps
}

((df - df.auto) / df)

df.auto <- matrix(0, dim(df)[1], dim(df)[2])
eps <- 1E-4

for(i in 1:length(pars)){
    pars.i1 <- pars.i2 <- pars
    pars.i1[i] <- pars[i] + eps
    pars.i2[i] <- pars[i] - eps
    Ires.i1 <- inneropt(data,times,pars.i1,coefs,lik,proc, in.meth,control.in)
    Ires.i2 <- inneropt(data,times,pars.i2,coefs,lik,proc, in.meth,control.in)
    ncoefs.i1 <- Ires.i1$coefs
    ncoefs.i2 <- Ires.i2$coefs
    devals.i1 = as.matrix(lik$bvals%*%ncoefs.i1)
    devals.i2 = as.matrix(lik$bvals%*%ncoefs.i2)
    colnames(devals.i1) <- colnames(devals.i2) <- proc$more$names
    f.i1 = as.vector( as.matrix(lik$more$fn(times, devals.i1, pars, lik$more$more))*sqrt(weights))
    f.i2 = as.vector( as.matrix(lik$more$fn(times, devals.i2, pars, lik$more$more))*sqrt(weights))
    df.auto[,i] <- (f.i1 - f.i2) / 2 / eps
}
((df - df.auto) / df)

## For d2fdcdp:
## H1
## H1.1:
H1.1a <- 2 * ((more$fn(more$qts, devals, pars, more$more) - ddevals) * 0 + pars["b"] * (1 - devals) * (- devals))
H1.1b <- 2 * ((more$fn(more$qts, devals, pars, more$more) - ddevals) * (1 - devals) + pars["b"] * (1 - devals)^2 * more$more$y.d)
cbind(H1.1[,1,1:2], H1.1a, H1.1b)
##H3:
H3m <- 2*((more$fn(more$qts, devals, pars, more$more)-ddevals) * pars["b"] * more$more$dy.d + (pars["b"] * more$more$y.d + pars["a"])* pars["b"] * (1 - devals) * more$more$dy.d)
##H3.1
H3.1m <- -2 * pars["b"] ^ 2 * (1 - devals) ^ 2 * more$more$dy.d
## H4m
H4m <- - 2 * pars["b"] * (1 - devals) * more$more$dy.d

##
eps <- 1E-3
delayProcObj1 <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, pars["tau"]+eps)
proc$more$more$y.d <- delayProcObj1$y.d
proc$more$more$dy.d <- delayProcObj1$dy.d
proc$more$more$bvals.d <- delayProcObj1$bvals.d
proc$more$more$dbvals.d <- delayProcObj1$dbvals.d
devals = as.matrix(proc$bvals$bvals %*% ncoefs)
ddevals = as.matrix(proc$bvals$dbvals %*% ncoefs)
dfdtau.m1 <- proc$more$fn(proc$more$qpts, devals, pars,proc$more$more)

as.vector(dfdtau.m1 - dfdtau.m2) / 2 / eps
## For d2fdc2
## H1.1:
H1.1m <- 2 * (more$fn(more$qts, devals, pars, more$more)-ddevals) * 0 + 2 * (pars["b"] * (1 - devals))^2
as.vector(H1.1[,1,] / H1.1m)

## H1.2
H1.2m <- 2 * (more$fn(more$qts, devals, pars, more$more)-ddevals) * (-pars["b"]) + 2 * pars["b"] * (1 - devals) * (- pars["b"] * more$more$y.d - pars["a"])
as.vector(H1.2[,1,] / H1.2m)

## Calculate \frac{\partial J}{\partial \psi} inside function SplineCoefsErr.DDE
dJdp <- dJdp.num <- matrix(0, length(proc$more$qpts), length(pars))
devals <- as.matrix(proc$bvals$bvals %*% coefs2)
ddevals <- as.matrix(proc$bvals$dbvals %*% coefs2)
colnames(devals) = proc$more$names
colnames(ddevals) = proc$more$names
nParsDelay <- sum(proc$more$more$tauIndex)

## Numerical:
eps <- 1E-4
for(i in 1:(length(pars) - nParsDelay)){
    pars.i1 <- pars.i2 <- pars
    pars.i1[i] <- pars[i] + eps
    pars.i2[i] <- pars[i] - eps
    J.i1 <- make.SSElik()$fn(ddevals, proc$more$qpts, devals, pars.i1, proc$more)
    J.i2 <- make.SSElik()$fn(ddevals, proc$more$qpts, devals, pars.i2, proc$more)
    dJdp.num[,i] <- (J.i1 - J.i2) / 2 / eps
}

for(i in (length(pars) - nParsDelay + 1):length(pars)){
    pars.i1 <- pars.i2 <- pars
    pars.i1[i] <- pars[i] + eps
    pars.i2[i] <- pars[i] - eps
    delayProcObj.i1 <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = pars.i1[i])
    proc$more$more$y.d <- delayProcObj.i1$y.d
    proc$more$more$dy.d <- delayProcObj.i1$dy.d
    proc$more$more$bvals.d <- delayProcObj.i1$bvals.d
    proc$more$more$dbvals.d <- delayProcObj.i1$dbvals.d
    J.i1 <- make.SSElik()$fn(ddevals, proc$more$qpts, devals, pars.i1, proc$more)
    delayProcObj.i2 <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = pars.i2[i])
    proc$more$more$y.d <- delayProcObj.i2$y.d
    proc$more$more$dy.d <- delayProcObj.i2$dy.d
    proc$more$more$bvals.d <- delayProcObj.i2$bvals.d
    proc$more$more$dbvals.d <- delayProcObj.i2$dbvals.d
    J.i2 <- make.SSElik()$fn(ddevals, proc$more$qpts, devals, pars.i2, proc$more)
    dJdp.num[,i] <- (J.i1 - J.i2)/ 2 / eps
}

## Analytical

proc$more$more$y.d <- delayProcObj$y.d
proc$more$more$dy.d <- delayProcObj$dy.d
proc$more$more$bvals.d <- delayProcObj$bvals.d
proc$more$more$dbvals.d <- delayProcObj$dbvals.d
fdevals = proc$more$fn(proc$more$qpts, devals, pars, proc$more$more)
difs <- ddevals - fdevals
weights <- checkweights(proc$more$weights, proc$more$whichobs, difs)

dfdp <- proc$more$dfdp(proc$more$qpts, devals, pars, proc$more$more)
dfdtau <-  proc$more$dfdtau(proc$qpts, devals, pars, proc$more)
dJdp[, 1:(length(pars) - nParsDelay)] <- sweep(dfdp[,1,], 1, 2 * difs[,1] * weights, FUN = "*")
dJdp[, (length(pars) - nParsDelay + 1):length(pars)] <- dfdtau[,1,] * 2 * difs[,1] * weights

dJdp.num / dJdp
