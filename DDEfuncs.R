ProfileSSE <- function(pars, allpars, times, data, coefs, lik, proc, in.meth = "nlminb",
                       control.in = NULL, active = 1:length(pars), dcdp = NULL,
                       oldpars = NULL, use.nls = TRUE, sgn = 1)
{
    allpars[active] = pars
    f <- ProfileSSE.AllPar(pars = allpars, times = times, data = data,
                           coefs = coefs, lik = lik, proc = proc, in.meth = in.meth,
                           control.in = control.in, dcdp = dcdp, oldpars = oldpars,
                           use.nls = use.nls, sgn = sgn)
    if (use.nls) {
    }
    else {
        f$gradient = f$gradient[, active, drop = FALSE]
        f$dcdp = f$dcdp[, active, drop = FALSE]
    }
    return(f)
}

Smooth.LS.DDE <- function(fn, data, times, pars, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", control.in = list(), eps = 1e-06, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL)
{
    dims = dim(data)
    ## if(dims[2] == 1){
    ##     data <- matrix(data[times > pars["tau"],], ncol = 1)
    ## }else
    ##     data <- data[times > pars["tau"],]
    ## coefs <- coefs[times > pars["tau"]]

    ## times <- times[times > pars["tau"]]

    profile.obj = LS.setup(pars, coefs, fn, basisvals, lambda,
        fd.obj, more, data, weights, times, quadrature, eps = 1e-06,
        posproc, poslik, discrete, names = names, sparse = sparse,
        likfn = make.id(), likmore = NULL)
    lik = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs

    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs0, "dimnames")[[2]]
    fdobj <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj, "class") <- "fd"
    y.d <- eval.fd(times - pars["tau"], fdobj)
    lik$more$y.d <- proc$more$y.d <- y.d

    Ires = inneropt.DDE(data, times, pars, coefs, lik, proc, in.meth,
        control.in)
    ncoefs = Ires$coefs
    Ires = Ires$res
    ncoefs = as.matrix(ncoefs)
    if (!is.null(proc$more$names)){
        colnames(ncoefs) = proc$more$names
    }
    if (!is.null(fd.obj)){
        if (length(dims) > 2) {
            ncoefs = array(ncoefs, c(length(ncoefs)/(dims[2] *
                dims[3]), dims[2], dims[3]))
        }
        else {
            ncoefs = array(ncoefs, c(length(ncoefs)/dims[2],
                dims[2]))
        }
        fd.obj = fd(ncoefs, fd.obj$basis)
        return(list(fd = fd.obj, lik = lik, proc = proc, inner.result = Ires))
    }
    else {
        return(list(coefs = ncoefs, lik = lik, proc = proc, inner.result = Ires,
            data = data, times = times))
    }
}

inneropt.DDE <- function(data, times, pars, coefs, lik, proc,
                         in.meth = "nlminb", control.in = list())
{
    check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    if (in.meth == "optim") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$maxit)) {
            control.in$maxit = 1000
        }
        if (is.null(control.in$reltol)) {
            control.in$reltol = 1e-12
        }
        if (is.null(control.in$meth)) {
            control.in$meth = "BFGS"
        }
        if (is.null(control.in$reportHessian)){
            control.in$reportHessian = TRUE
        }
        imeth = control.in$meth
        control.in$meth = NULL
        res = optim(coefs, SplineCoefsErr, gr = SplineCoefsDC,
            hessian = control.in$reportHessian, control = control.in,
            times = times, data = data, lik = lik, proc = proc,
            pars = pars, y.d = y.d, method = imeth)
        ncoefs = matrix(res$par, ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
    }
    else if (in.meth == "nlminb") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$eval.max)) {
            control.in$eval.max = 2000
        }
        if (is.null(control.in$iter.max)) {
            control.in$iter.max = 1000
        }
        if (is.null(control.in$rel.tol)) {
            control.in$rel.tol = 1e-12
        }
        if (is.null(control.in$useHessian)) {
            Hessian = SplineCoefsDC2
        }
        else {
            Hessian = NULL
        }
        ## SplineCoefsErr do not need to be changed.
        ##
        res = nlminb(coefs, SplineCoefsErr, gradient = SplineCoefsDC,
            hessian = Hessian, control = control.in, times = times,
            data = data, lik = lik, proc = proc, pars = pars, y.d = y.d)
        ncoefs = matrix(res$par, ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
    }
    else {
        stop("Unknown optimizer specified")
    }
    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
    return(list(coefs = ncoefs, res = res))
}

ProfileSSE.AllPar.DDE <- function(pars, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1,
                                  tauIndex = NULL, basisObj)
{
    ## Squared Error outer criterion

    ##    coefs = as.vector(coefs)  # First run the inner optimization
    ##
    f1 = SplineCoefsErr(coefs,times,data,lik,proc,pars)
    tau <- pars[tauIndex]
    if(use.nls){
        # If using NLS, we need to keep track
        if(file.exists('curcoefs.tmp')){
          altcoefs = as.matrix(read.table('curcoefs.tmp'))
          if( !(length(altcoefs)==length(coefs)) ){
              stop(paste('Variables in curcoefs.tmp do not conform;',
                       'file exists from previous experiments?'))
          }
      } else {
          altcoefs = coefs
        }
        if(file.exists('counter.tmp')){
          counter = read.table('counter.tmp')
          niter = counter[nrow(counter),1]
        } else {
          counter = matrix(c(1,0,pars),1,length(pars)+2)
          niter = 0
        }

        f2 = SplineCoefsErr(altcoefs,times,data,lik,proc,pars)

        if(f2 < f1){
            coefs = altcoefs
            f1 = f2
        }

        altdevals = as.matrix(lik$bvals %*% matrix(altcoefs,ncol(lik$bvals),
                              length(altcoefs)/ncol(lik$bvals)))
        colnames(altdevals) = proc$more$names
    }

    ## The coefficients for delayed times:
    bvals.d <- eval.basis(times-tau, basisObj)
    y.d <- eval.fd(times - tau, )
    ##################################################
    ## Not Sure:
    ## if(!is.null(dcdp)){
    ##     tcoefs = as.vector(coefs) + dcdp%*%(pars-oldpars);
    ##     f2 = SplineCoefsErr(tcoefs,times,data,lik,proc,pars)
    ##     if(f2 < f1){
    ##         coefs = tcoefs
    ##         f1 = f2
    ##     }
    ## }
    ###################################################

    ## Inner optimization need to be changed as well.
    Ires = inneropt(data,times,pars,coefs,lik,proc,y.d, in.meth,control.in)
    ncoefs = Ires$coefs

    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names

    ## Squared errors:
    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
    f = as.vector( as.matrix(data -
        lik$more$fn(times, devals, pars, lik$more$more))*sqrt(weights))
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
    ## Changed for DDE: xd stands for delayed x.
    d2Hdc2  = SplineCoefsDC2sparse.DDE(ncoefs,times,data,lik,proc,pars, bvalsd = bvalsd)
    d2Hdcdp = SplineCoefsDCDP(ncoefs,times,data,lik,proc,pars)
    ## Got warning message:
    ## In dim(weights[whichrows, ]) == dim(diffs) :
    ## longer object length is not a multiple of shorter object length

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

    if(!is.null(lik$report)){ print(f) }

    f = sgn*f
    df = sgn*df

    if(use.nls){
        tf = sum(lik$fn(data,times,devals,pars,lik$more))
        tf2 = sum(lik$fn(data,times,altdevals,pars,lik$more))

        if(tf <= tf2){
            write.table(ncoefs,file='curcoefs.tmp',
                        row.names=FALSE,col.names=FALSE)
            niter = counter[nrow(counter),1]
        }

        if(niter==0){
          counter[1,2] = tf
          write.table(counter,file='counter.tmp',
                      col.names=FALSE,row.names=FALSE)
        }

     if(niter > 1){
      if(tf < counter[niter,2]){
       counter = rbind(counter,c(niter+1,tf,pars))
       write.table(counter,file='counter.tmp',col.names=FALSE,row.names=FALSE)
      }
    }

        write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
        attr(f,'gradient') = df
        return(f)
    }
    else{
        return(list(value=f,gradient=df,coefs=ncoefs,dcdp=dcdp))
    }
}

## Do not need to use this function. Unchanged for DDE.
SplineCoefsErr.DDE <- function(coefs, times, data, lik, proc, pars, sgn = 1){
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2,proc$bvals, pars, proc$more)
    if (!is.null(proc$report)){
        print(f)
    }
    return(sgn * f)
}

SplineCoefsDC.DDE <- function (coefs, times, data, lik, proc, pars, sgn = 1)
{
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    ## lik$dfdx is not changed
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals,
        pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars,
        proc$more)
    g = as.vector(g)
    return(sgn * g)
}


## Needs changes
proc$delay$dfdc <- function(coefs, bvals, pars, more, delay){
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    g1 <- delay$dfdx(ddevals, more$qpts, devals, pars, more)
    weights = checkweights(more$weights, more$whichobs, g1)
    g2 = weights * (ddevals - more$fn(more$qpts, devals, pars, more$more))
    g = as.vector(as.matrix(t(bvals$bvals) %*% g1 + 2 * t(bvals$dbvals) %*% g2))
    return(g)

}

SplineCoefsDC2sparse.DDE <- function(coefs,times,data,lik,proc,pars,sgn=1, bvalsd)
{
    ## This function will not be changed.
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

## "more" is taken from $proc$more
## Originally, proc$d2fdc2() has the result of
## $$ \bigl[\dot{x}_{l}(t)-f_{l}(t)\bigr]\bigl[\frac{df_{l}}{dx_{j}dx_{i}}
## + \frac{df_{l}}{dx_{j}} \bigl[\frac{df_{l}}{dx_{i}}\bigr]^{\top} $$

d2fdc2.DDE <- function(coefs, bvals, pars, more){
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    ## H1: make.SSElik()$d2fdx2 is the same function as lik$d2fdx2 ?!
    ## But the arguments may be different depending on the
    ## Here more is proc$more
    H1 = delay$d2fdx2(ddevals, more$qpts, devals, pars,
        more, devals.d, ddevals.d)
    H2 = delay$dfdx(more$qpts, devals, pars, more$more, devals.d)
    weights = checkweights(more$weights, more$whichobs, H1[,
        , 1, drop = TRUE])

    ##################################################
    ## New for DDE:
    ##################################################
    H2.1 <- delay$dfdx.d(more$qpts, devals, pars, more$more, devals.d)
    H1.1 <- delay$d2fdx.d2(ddevals, more$qpts, devals, pars, more, devals.d, ddevals.d)
    H1.2 <- delay$d2fdxdx.d(ddevals, more$qpts, devals, pars, more, devals.d, ddevals.d)
    H1.3 <- delay$d2fdx.ddx(ddevals, more$qpts, devals, pars, more, devals.d, ddevals.d)
    H = list(len = dim(bvals$bvals)[2])
    for (i in 1:dim(devals)[2]){
        H[[i]] = list(len = dim(devals))
        for (j in 1:dim(devals)[2]) {
            H[[i]][[j]] <- t(bvals$bvals) %*% diag(H1[, i, j]) %*% bvals$bvals +
                t(bvals$bvals.d) %*% diag(H1.1[, i, j]) %*% bvals$bvals.d +
                t(bvals$dbvals) %*% diag(H1.2[, i, j]) %*% bvals$bvals +
                t(bvals$bvals) %*% diag(H1.2[, j, i]) %*% bvals$dbvals -
                2 * t(bvals$dbvals) %*% diag(H2[, i, j] * weights[, i]) %*% bvals$bvals -
                2 * t(bvals$bvals) %*% diag(H2[, j, i] * weights[, j]) %*% bvals$dbvals -
                2 * t(bvals$dbvals) %*% diag(H2.1[, i, j] * weights[, i]) %*% bvals$bvals.d -
                2 * t(bvals$bvals.d) %*% diag(H2.1[, j, i] * weights[, j]) %*% bvals$dbvals
        }
        H[[i]][[i]] = H[[i]][[i]] + 2 * t(bvals$dbvals) %*% diag(weights[,
            i]) %*% bvals$dbvals
    }
    H = blocks2mat(H)
    return(H)
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
proc$d2fdcdp.d <- function (coefs, bvals, pars, more)
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
    ##################################################
    ## New for DDE:
    ##################################################
    H1.1 <- delay$d2fdxdp(ddevals, more$qpts, devals, pars, more)
    H3 <- delay$d2fdxdd(ddevals, more$qpts, devals, pars, more)
    H3.1 <- delay$d2fdxdd.d(ddevals, more$qpts, devals, pars, more, devals.d)
    H4 <- 2 * delay$dfdd(more$qpts, devals, pars, more$more, devals.d)
    H = c()
    for (i in 1:(length(pars)-nParsDelay)) {
        H = cbind(H, as.vector(as.matrix(t(bvals$bvals) %*% H1[,, i] + t(bvals$bvals.d) %*% H1.1[,,i] - t(bvals$dbvals) %*% (weights * H2[, , i]))))
    }
    for(i in (length(pars) + 1 - nParsDelay): length(pars)){
        H = cbind(H, as.vector(as.matrix(t(bvals$bvals) %*% H3[,, i] + t(bvals$bvals.d) %*% H3.1[,,i] - t(bvals$dbvals) %*% (weights * H4[, , i]))))
    }
    return(H)
}

delay$d2fdx2(ddevals, more$qpts, devals, pars, more, devals.d, ddevals.d)

## Replacing make.SSElik()$dfdx()
proc$delay$dfdx <- function(data, times, devals, pars, more){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdx.d = more$delay$dfdx.d(times, devals, pars, more$more, y.d)
    g = c()
    for (i in 1:dim(dfdx)[3]) {
        g = cbind(g, apply(difs * (dfdx[, , i] + dfdx.d[, , i]), 1, sum))
    }
    return(-2 * g)
}


weights = checkweights(more$weights, more$whichobs, H1[,, 1, drop = TRUE])

delay$d2fdx.d2 <- function(){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx.d = delay$dfdx.d(times, devals, pars, more$more)
    d2fdx.d2 = delay$d2fdx.d2(times, devals, pars, more$more, devals.d)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdx2.d)[3]) {
        for (j in 1:dim(d2fdx2.d)[4]) {
            H[, i, j] = apply(-difs * d2fdx2.d[, , i, j] + weights *
             dfdx.d[, , j] * dfdx.d[, , i], 1, sum)
        }
    }
    return(2 * H)
}

delay$d2fdxdx.d <- function(){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdx.d <- more$dfdx(times, devals, pars, more$more)
    d2fdx2 <- delay$d2fdxdx.d(times, devals, pars, more$more, devals.d)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdx2)[3]) {
        for (j in 1:dim(d2fdx2)[4]) {
            H[, i, j] = apply(-difs * d2fdx2[, , i, j] + weights *
             dfdx[, , j] * dfdx.d[, , i], 1, sum)
        }
    }
    return(2 * H)
}

delay$d2fdx.ddx <- function(){
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdx.d <- more$dfdx(times, devals, pars, more$more)
    d2fdx2 <- delay$d2fdx.ddx(times, devals, pars, more$more, devals.d)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, difs)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdx2)[3]){
        for (j in 1:dim(d2fdx2)[4]){
            H[, i, j] = apply(-difs * d2fdx2[, , i, j] + weights *
             dfdx.d[, , j] * dfdx[, , i], 1, sum)
        }
    }
    return(2 * H)
}

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
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
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
vectorD <- list()

## y.d is added for delay
vectorD.Fun <- vectorD$proc$more$fn <- function (times, y, p, more)
{
    y.d <- more$y.d
    r <- y
    r[, "S"] <- p["b"] * y.d[,"S"] * (1 - y[,"S"]) - p["a"] * y[,"S"]
    return(r)
}

vectorD$proc$more$dfdx <- function(t, x, pars, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[,"S", "S"] <- -p["b"] * y.d[,"S"] - p["a"]
    return(r)
}

vectorD$proc$more$dfdp <- function(times, y, p, more)
{
    y.d <- more$y.d
    r = array(0, c(dim(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    r[, "V", "a"] <- - y[,"S"]
    r[, "V", "b"] <- y.d[,"S"] * (1 - y[,"S"])
    return(r)
}

## Return an array of all 0's
vectorD$proc$more$d2fdx2 <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

vectorD$proc$more$d2fdxdp <- function (times, y, p, more)
{
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2], length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[,"S", "S", "a"] <- -1
    r[,"S", "S", "b"] <- - y.d[,"S"]
}

## New functions:
vectorD$proc$more$dfdx.d <- function(times, x, pars, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y.d))
    r[,"S", "Delay"] <- p["b"] * (1 - y[,"S"])
    return(r)
}

vectorD$proc$more$d2fdxdx.d <- function(times, x, pars, more){
    y.d <- more$y.d
    r = array(-1, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y.d))
    return(r)
}

vectorD$proc$more$d2fdxdx.d <- function(times, x, pars, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y.d), colnames(y.d))
    return(r)
}

vectorFun <- list()

vectorFun$fn <- function (times, y, p, more)
{
    y.d <- more$y.d
    r <- y
    dimnames(r) <- dimnames(y) ## ??
    r[, "S"] <- p["b"] * y.d[,"S"] * (1 - y[,"S"]) - p["a"] * y[,"S"]
    return(r)
}

vectorFun$fn.ode <- function(times, y, p)
{
    r <- y
    dimnames(r) <- dimnames(y) ## ??
    r[, "S"] <- p["b"] * y.d[,"S"] * (1 - y[,"S"]) - p["a"] * y[,"S"]
    return(r)
}


vectorFun$dfdx <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[,"S", "S"] <- -p["b"] * y.d[,"S"] - p["a"]
    return(r)
}

vectorFun$dfdp <- function(times, y, p, more)
{
    y.d <- more$y.d
    r = array(0, c(dim(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    r[, "V", "a"] <- - y[,"S"]
    r[, "V", "b"] <- y.d[,"S"] * (1 - y[,"S"])
    return(r)
}

vectorFun$d2fdxdp <- function (times, y, p, more)
{
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2], length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[,"S", "S", "a"] <- -1
    r[,"S", "S", "b"] <- - y.d[,"S"]
}

vectorFun$d2fdx2 <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

vectorFun$dfdx.d <- function(times, x, pars, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y.d))
    r[,"S", "Delay"] <- p["b"] * (1 - y[,"S"])
    return(r)
}

vectorFun$d2fdxdx.d <- function(times, x, pars, more){
    y.d <- more$y.d
    r = array(-1, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y.d))
    return(r)
}

vectorFun$d2fdxdx.d <- function(times, x, pars, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y.d), colnames(y.d))
    return(r)
}
