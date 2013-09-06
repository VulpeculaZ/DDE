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
