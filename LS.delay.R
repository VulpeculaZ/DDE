LS.setup.DDE <- function(pars, coefs=NULL, fn, basisvals=NULL, lambda, fd.obj=NULL,
                     more=NULL, data=NULL, weights=NULL, times=NULL,
                     quadrature=NULL, eps=1e-6, posproc=FALSE, poslik = FALSE,
                     discrete=FALSE, names=NULL, sparse=FALSE,
                     likfn = make.id(), likmore = NUL)
{
    ## Functions replaced:
    ## proc$dfdc need to be replaced by the new proc$dfdc
    LS.result <- LS.setup(pars, coefs, fn, basisvals, lambda, fd.obj,
                     more, data, weights, times,
                     quadrature, eps, posproc, poslik,
                     discrete, names, sparse,
                     likfn = make.id(), likmore)

}

delay <- list()

## Needs changes
## Replace the original proc$dfdc
delay$dfdc <- function(coefs, bvals, pars, more){
    ## delay <- more$delay
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    ## Need change:
    ## dfdx.d need to be tested.
    g1 <- make.SSElik()$dfdx(ddevals, more$qpts, devals, pars, more)
    g1.d <- more$delay$dfdx.d(ddevals, more$qpts, devals, pars, more)
    weights = checkweights(more$weights, more$whichobs, g1)
    g2 = weights * (ddevals - more$fn(more$qpts, devals, pars, more$more))
    bvals.d <- rbind( matrix(0, dim(g1.d)[1] - dim(more$more$bvals.d)[1], dim(g1.d)[1]), more$more$bvals.d)
    g = as.vector(as.matrix(t(bvals$bvals) %*% g1 + t(bvals.d) %*% g1.d + 2 * t(bvals$dbvals) %*% g2))
    return(g)
}

## Replacing make.SSElik()$dfdx() for calculating g1.d in the new proc$dfdc
## put into proc$more$delay$dfdx.d
delay$dfdx.d <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx.d = more$dfdx.d(times, devals, pars, more$more)
    g = c()
    for (i in 1:dim(dfdx.d)[3]) {
        g = cbind(g, apply(difs * dfdx.d[, , i], 1, sum))
    }
    return(-2 * g)
}


## "more" is taken from $proc$more
## To replace the original proc$d2fdc2() has the result of
## $$ \bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\bigl[\frac{df_{l}}{dx_{j}dx_{i}}
## + \frac{df_{l}}{dx_{j}} \bigl[\frac{df_{l}}{dx_{i}}\bigr]^{\top} $$

delay$d2fdc2.DDE <- function(coefs, bvals, pars, more){
    delay <- more$delay
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    ## H1: make.SSElik()$d2fdx2 is the same function as lik$d2fdx2 ?!
    ## But the arguments may be different depending on the
    ## Here more is proc$more
    H1 = make.SSElik()$d2fdx2(ddevals, more$qpts, devals, pars,
        more)
    H2 = more$dfdx(more$qpts, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, H1[,
        , 1, drop = TRUE])
    ##################################################
    ## New for DDE:
    ##################################################
    H2.1 <- more$dfdx.d(more$qpts, devals, pars, more$more)
    H1.1 <- delay$d2fdx.d2(ddevals, more$qpts, devals, pars, more)
    H1.2 <- delay$d2fdxdx.d(ddevals, more$qpts, devals, pars, more)
    ## H1.3 <- delay$d2fdx.ddx(ddevals, more$qpts, devals, pars, more)
    bvals.d <- rbind( matrix(0, length(H1[,1,1]) - dim(more$more$bvals.d)[1], length(H1[,1,1])), more$more$bvals.d)
    dbvals.d <- rbind( matrix(0, length(H1[,1,1]) - dim(more$more$dbvals.d)[1], length(H1[,1,1])), more$more$dbvals.d)
    H = list(len = dim(more$more$bvals)[2])
    for (i in 1:dim(devals)[2]){
        H[[i]] = list(len = dim(devals))
        for (j in 1:dim(devals)[2]){
            H[[i]][[j]] <- t(bvals$bvals) %*% diag(H1[, i, j]) %*% bvals$bvals +
                t(bvals.d) %*% diag(H1.1[, i, j]) %*% bvals.d +
                t(bvals$bvals) %*% diag(H1.2[, i, j]) %*% bvals.d +
                t(bvals.d) %*% diag(H1.2[, j, i]) %*% bvals$bvals -
                2 * t(bvals$dbvals) %*% diag(H2[, i, j] * weights[, i]) %*% bvals$bvals -
                2 * t(bvals$bvals) %*% diag(H2[, j, i] * weights[, j]) %*% bvals$dbvals -
                2 * t(bvals$dbvals) %*% diag(H2.1[, i, j] * weights[, i]) %*% bvals.d -
                2 * t(bvals.d) %*% diag(H2.1[, j, i] * weights[, j]) %*% bvals$dbvals
        }
        H[[i]][[i]] = H[[i]][[i]] + 2 * t(bvals$dbvals) %*% diag(weights[,
            i]) %*% bvals$dbvals
    }
    H = blocks2mat(H)
    return(H)
}

##################################################
## Used in proc$d2fdc2:
##    H2.1 <- more$dfdx.d
##    H1.1 <- delay$d2fdx.d2
##    H1.2 <- delay$d2fdxdx.d
##    H1.3 <- delay$d2fdx.ddx
##################################################

delay$d2fdx.d2 <- function(data, times, devals, pars, more ){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx.d = more$dfdx.d(times, devals, pars, more$more)
    d2fdx.d2 = more$d2fdx.d2(times, devals, pars, more$more)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdx.d2)[3]) {
        for (j in 1:dim(d2fdx.d2)[4]) {
            H[, i, j] = apply(-difs * d2fdx.d2[, , i, j] + weights *
             dfdx.d[, , i] * dfdx.d[, , j], 1, sum)
        }
    }
    return(2 * H)
}

delay$d2fdxdx.d <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdx.d <- more$dfdx.d(times, devals, pars, more$more)
    d2fdxdx.d <- more$d2fdxdx.d(times, devals, pars, more$more)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdxdx.d)[3]) {
        for (j in 1:dim(d2fdxdx.d)[4]) {
            H[, i, j] = apply(-difs * d2fdxdx.d[, , i, j] + weights * (dfdx[, , i] * dfdx.d[, , j]), 1, sum)
        }
    }
    return(2 * H)
}

delay$d2fdx.ddx <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdx.d <- more$dfdx.d(times, devals, pars, more$more)
    d2fdxdx.d <- more$d2fdxdx.d(times, devals, pars, more$more)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdx2)[3]){
        for (j in 1:dim(d2fdx2)[4]){
            H[, i, j] = apply(-difs * d2fdxdx.d[, , i, j] + weights *
             dfdx.d[, , i] * dfdx[, , j], 1, sum)
        }
    }
    return(2 * H)
}


## To replace the original proc$d2fdcdp
## need to be debugged.
delay$d2fdcdp.DDE <- function (coefs, bvals, pars, more)
{
    nParsDelay <- sum(more$more$tauIndex)
    delay <- more$delay
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    bvals.d <- rbind(matrix(0, nrow = dim(bvals$bvals)[1] - dim(more$more$bvals.d)[1], ncol = dim(bvals$bvals)[2]), more$more$bvals.d)
    dbvals.d <- rbind(matrix(0, nrow = dim(bvals$dbvals)[1] - dim(more$more$dbvals.d)[1], ncol = dim(bvals$dbvals)[2]), more$more$dbvals.d)
    ## Should I use the full length dy.d?
    ## ddevals.d <- as.matrix(more$more$dbvals.d %*% coefs)
    ## ddevals.d <- more$more$dy.d
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    H1 = make.SSElik()$d2fdxdp(ddevals, more$qpts, devals, pars,
        more)
    H2 = 2 * more$dfdp(more$qpts, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, H1[,
        , 1, drop = FALSE])
    ##################################################
    ## New for DDE:
    ##################################################
    H1.1 <- delay$d2fdx.ddp(ddevals, more$qpts, devals, pars, more)
    H3 <- delay$d2fdxdtau(ddevals, more$qpts, devals, pars, more)
    H3.1 <- delay$d2fdx.ddtau(ddevals, more$qpts, devals, pars, more)
    H3.2 <- delay$d2fdx.ddtau2(ddevals, more$qpts, devals, pars, more)
    H4 <- 2 * more$dfdtau(more$qpts, devals, pars, more)[,,more$more$tauIndex, drop = FALSE]
    H = c()
    for (i in 1:(length(pars)-nParsDelay)) {
        H = cbind(H, as.vector(as.matrix(t(bvals$bvals) %*% H1[,, i] + t(bvals.d) %*% H1.1[,,i] - t(bvals$dbvals) %*% (weights * H2[, , i]))))
    }
    for(i in 1: nParsDelay){
        H <- cbind(H, as.vector(as.matrix(t(bvals$bvals) %*% H3[,, i] + t(bvals.d) %*% H3.1[,,i] - t(dbvals.d) %*% H3.2[,,i]- t(bvals$dbvals) %*% (weights * H4[, , i]))))
    }
    return(H)
}

##################################################
## Funtions to be used in proc$d2fdcdp.DDE
## H1.1 <- delay$d2fdx.ddp
## H3 <- delay$d2fdxdtau
## H3.1 <- delay$d2fdx.ddtau
## H4 <- 2 * more$dfdtau
##################################################
## More is from proc$more$delay



## H1.1 <- delay$d2fdx.ddp() :
## -2\bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\frac{d^{2}f_{l}}{dx_{j}d\theta^{\top}}+2\frac{df_{l}}{dx_{j}}\frac{df_{l}}{d\theta^{\top}}
delay$d2fdx.ddp <- function (data, times, devals, pars, more)
{
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx.d = more$dfdx.d(times, devals, pars, more$more)
    dfdp = more$dfdp(times, devals, pars, more$more)
    d2fdx.ddp = more$d2fdx.ddp(times, devals, pars, more$more)
    H = array(0, c(dim(devals), length(pars)-sum(more$tauIndex)))
    for (i in 1:dim(d2fdx.ddp)[3]) {
        for (j in 1:dim(d2fdx.ddp)[4]) {
            H[, i, j] = apply(-difs * d2fdx.ddp[, , i, j] + weights *
                dfdx.d[, , i] * dfdp[, , j], 1, sum)
        }
    }
    return(2 * H)
}

## H3 <- delay$d2fdxdtau() :
## 2\frac{df_{l}}{dx_{j}}\frac{df_{l}}{d\tau^{\top}}-2\bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\frac{d^{2}f_{l}}{dx_{j}d\tau^{\top}}
delay$d2fdxdtau <- function (data, times, devals, pars, more)
{
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdtau <- more$dfdtau(times, devals, pars, more)[,,more$more$tauIndex, drop = FALSE]
    d2fdxdtau <- more$d2fdxdtau(times, devals, pars, more)[,,,more$more$tauIndex, drop = FALSE]
    H = array(0, c(dim(devals), sum(more$more$tauIndex)))
    for (i in 1:dim(d2fdxdtau)[3]) {
        for (j in 1:dim(d2fdxdtau)[4]) {
            H[, i, j] = apply(-difs * d2fdxdtau[, , i, j] + weights *
                dfdx[, , i] * dfdtau[, , j], 1, sum)
        }
    }

    return(2 * H)
}

## H3.1 <- delay$d2fdx.ddtau() :
## 2\frac{df_{l}}{dx_{j}(-T_{j})}\frac{df_{l}}{d\tau^{\top}}-2\bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\frac{d^{2}f_{l}}{dx_{j}(-T_{j})d\tau^{\top}}
delay$d2fdx.ddtau <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx.d = more$dfdx.d(times, devals, pars, more$more)
    dfdtau <- more$dfdtau(times, devals, pars, more)[,,more$more$tauIndex, drop = FALSE]
    d2fdx.ddtau <- more$d2fdx.ddtau(times, devals, pars, more)[,,, more$more$tauIndex, drop = FALSE]
    H = array(0, c(dim(devals), sum(more$more$tauIndex)))
    for (i in 1:dim(d2fdx.ddtau)[3]) {
        for (j in 1:dim(d2fdx.ddtau)[4]) {
            H[, i, j] = apply(-difs * d2fdx.ddtau[, , i, j] + weights *
                dfdx.d[, , i] * dfdtau[, , j], 1, sum)
        }
    }
    return(2 * H)
}

## H3.2
delay$d2fdx.ddtau2 <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx.d <- more$dfdx.d(times, devals, pars, more$more)[,,more$more$tauIndex, drop = FALSE]
    H = array(0, c(dim(devals), sum(more$more$tauIndex)))
    for (i in 1:dim(H)[3]) {
            H[, i, i] = apply(-difs * weights *
                dfdx.d[, , i], 1, sum)
    }
    return(2 * H)
}

##################################################
## more contains function to calculate derivatives: dfdx, dfdx.d and spline coefs.
##################################################

## H4 <- 2 * more$dfdtau() :
## 2\frac{df_{l}}{d\tau^{\top}}
dfdtau.DDE <- function(times, y, p, more){
    dx.d <- more$more$dy.d
    dfdx.d <- more$dfdx.d(times, y, p, more$more)
    r <- array(0, dim = dim(dfdx.d))
    for(i in 1:dim(dfdx.d)[2]){
        r[,i,] <-  -dfdx.d[,i,] * dx.d
    }
    dimnames(r) <- dimnames(dfdx.d)
    return(r)
}

d2fxdtau.DDE <- function(times, y, p, more){
    dx.d <- more$more$dy.d
    d2fdxdx.d <- more$d2fdxdx.d(times, y, p, more$more)
    r <- dx.dArr <- array(0, dim = dim(d2fdxdx.d))
    for(i in 1:dim(r)[2]){
        for(j in 1:dim(r)[3]){
            r[,i,j,] <- - d2fdxdx.d[,i,j,] * dx.d
        }
    }
    dimnames(r) <- dimnames(d2fdxdx.d)
    return(r)
}

d2fdx.ddtau.DDE <- function(times, y, p, more){
    dx.d <- more$more$dy.d
    d2fdx.d2 <- more$d2fdx.d2(times, y, p, more$more)
    r <- array(0, dim = dim(d2fdx.d2))
    for(i in 1:dim(r)[2]){
        for(j in 1:dim(r)[3]){
            r[,i,j,] <- - d2fdx.d2[,i,j,] * dx.d
        }
    }
    dimnames(r) <- dimnames(d2fdx.d2)
    return(r)
}
