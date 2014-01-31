## Need to use changed proc$dfdc
## No changes for this one!
## SplineCoefsDC.DDE <- function (coefs, times, data, lik, proc, pars, sgn = 1)
## {
##     coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
##     devals = as.matrix(lik$bvals %*% coefs2)
##     colnames(devals) = proc$more$names
##     ## lik$dfdx is not changed
##     g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals,
##         pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars,
##         proc$more)
##     g = as.vector(g)
##     return(sgn * g)
## }


## Have not been changed
## Use the new proc$d2fdc2 to replace the old one.
SplineCoefsDC2sparse.DDE <- function(coefs,times,data,lik,proc,pars,sgn=1)
{

    ## This function will not be changed after this point.
    ##  Inner Hessian with respect to coefficients
    coefs2 = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs2)
    colnames(devals) = proc$more$names

    ## H[,i,i] will be identity matrix
    d2lik = lik$d2fdx2(data,times,devals,pars,lik$more)

    H = list(len=ncol(lik$bvals))
    for(i in 1:dim(d2lik)[2]){
      H[[i]] = list(len=ncol(devals))
        for(j in 1:dim(d2lik)[3]){
            H[[i]][[j]] <- t(lik$bvals)%*%diag(d2lik[,i,j])%*%lik$bvals
        }
    }
    ## This will not be changed for DDE
    H = blocks2mat(H)
    H = H + proc$d2fdc2(coefs2,proc$bvals,pars,proc$more)
    return(sgn*H)
}

## Use the delay$d2fdcdp.DDE to replace the original proc$d2fdcdp
SplineCoefsDCDP.DDE <- function (coefs, times, data, lik, proc, pars, sgn = 1)
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
    ## Till now, H has all 0's
    H = H + proc$d2fdcdp(coefs2, proc$bvals, pars, proc$more)
    return(as.matrix(sgn * H))
}

##

d2fdxdp.tmp <- function (data, times, devals, pars, more)
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
