
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Evaluate fd objects at delayed times.
##' @param fd0 fd object starting at time 0 and ends at time tau
##' @param fd.d fd object starting at time tau
##' @param times Time points at which solution is evaluated. The values of the solution at times - tau are returned.
##' @param tau The delay
##' @param tau.max The Maximum of the delay.
##' @return A vector of values of the solution at times - tau are returned.
##' @author Ziqian Zhou
delay.fit <- function(fd0, fd.d, times, tau){
    basisvals0 <- fd0$basis
    basisvals.d <- fd.d$basis
    times.d <- times - tau
    start.d <- fd.d$basis$rangeval[1]
    y.d0 <- eval.fd(times.d[times.d < start.d], fd0) ## Same as bvals %*$ coefs
    dy.d0 <- eval.fd(times.d[times.d < start.d], fd0, Lfdobj = 1)
    y.d.d <- eval.fd(times.d[times.d >= start.d], fd.d)
    dy.d.d <- eval.fd(times.d[times.d >= start.d], fd.d, Lfdobj = 1)
    y.d <- rbind(y.d0, y.d.d)
    dy.d <- rbind(dy.d0, dy.d.d)
    ## Need to set \frac{dx(t_{i}-T_{j})}{d\mathbf{c}}=0 for t_{i}<T_{j}
    ## !!
    bvals.d <- eval.basis(times.d[times.d >= start.d], basisvals.d, 0)
    dbvals.d <- eval.basis(times.d[times.d >= start.d], basisvals.d, 1)
    return(list(y.d = y.d, dy.d = dy.d, dbvals.d = dbvals.d, bvals.d = bvals.d))
}


ProfileSSE.DDE <- function(pars, allpars, times, data, coefs, lik, proc, in.meth = "nlminb", control.in = NULL, active = 1:length(pars), dcdp = NULL, oldpars = NULL, use.nls = TRUE, sgn = 1, basisvals = NULL, fdobj0 = NULL)
{
    allpars[active] = pars
    f <- ProfileSSE.AllPar.DDE(pars = allpars, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth,control.in = control.in, dcdp = dcdp, oldpars = oldpars, use.nls = use.nls, sgn = sgn, basisvals = basisvals, fdobj0 = fdobj0, )
    if (use.nls) {
        attr(f,'gradient') = attr(f, 'gradient')[,active]
    }
    else{
        f$gradient = f$gradient[, active, drop = FALSE]
        f$dcdp = f$dcdp[, active, drop = FALSE]
    }
    return(f)
}

Profile.LS.DDE <- function(fn, data, times, pars, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, delay = NULL, tauMax = NULL,
    basisvals0 = NULL, coefs0 = NULL, tauIndex = NULL)
{
    if (is.null(active)) {
        active = 1:length(pars)
    }
    ##################################################
    ## Added to orginal Profile.LS
    if(is.null(tauMax)) tauMax <- min(times) + 1/3 * (range(times)[2]- range(times)[1])
    ## Prepare to calculate the delay:
    ## data.d <- data[times>=tauMax, ,drop = FALSE]
    ## times.d <- knots.d <- times[times >= tauMax]
    ## norder <- basisvals$nbasis - length(basisvals$params)
    ## nbasis.d <- length(knots.d) + norder - 2
    ## range.d <- range(knots.d)
    ## basisvals.d and basisvals.0 are supplied rather than created
    ## basisvals.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)
    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"

    profile.obj = LS.setup(pars = pars, coefs = coefs, fn = fn,
        basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times

    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = pars["tau"])
    delayLikObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = times, tau = pars["tau"])
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$dy.d <- delayLikObj$dy.d
    proc$more$more$dy.d <- delayProcObj$dy.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    lik$more$more$dbvals.d <- delayLikObj$dbvals.d
    proc$more$more$dbvals.d <- delayProcObj$dbvals.d
    proc$more$more$tauIndex <- tauIndex

    proc$dfdc <- delay$dfdc
    proc$d2fdc2 <- delay$d2fdc2.DDE
    proc$d2fdcdp <- delay$d2fdcdp.DDE
    proc$more$delay <- delay
    proc$more$dfdtau <- dfdtau.DDE
    proc$more$d2fdxdtau <- d2fxdtau.DDE
    proc$more$d2fdx.ddtau <- d2fdx.ddtau.DDE

    if (file.exists("curcoefs.tmp")) {
        file.remove("curcoefs.tmp")
    }
    if (file.exists("optcoefs.tmp")) {
        file.remove("optcoefs.tmp")
    }
    if (file.exists("counter.tmp")) {
        file.remove("counter.tmp")
    }
    Ires <- inneropt.DDE(data, times, pars, coefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ## Ires <- IresTmp
    ncoefs = Ires$coefs
    write.table(ncoefs, file = "optcoefs.tmp", col.names = FALSE,
        row.names = FALSE)
    write.table(ncoefs, file = "curcoefs.tmp", col.names = FALSE,
        row.names = FALSE)
    apars = pars[active]
    aparamnames = names(apars)
    if (out.meth == "ProfileGN") {
        res = Profile.GausNewt(pars = pars, times = times, data = data,
            coefs = ncoefs, lik = lik, proc = proc, in.meth = in.meth,
            control.in = control.in, active = active, control = control.out)
        apars = res$pars[active]
        ncoefs = res$in.res$coefs
        g = res$in.res$df
        resid = res$in.res$f
    }
    if (out.meth == "nls") {
        if (is.null(control.out$trace)) {
            control.out$trace = TRUE
        }
        if (is.null(control.out$maxiter)) {
            control.out$maxiter = 100
        }
        if (is.null(control.out$tol)) {
            control.out$tol = 1e-08
        }
        if (is.null(control.out$printEval)) {
            control.out$printEval = TRUE
        }
        if (is.null(control.out$warnOnly)) {
            control.out$warnOnly = TRUE
        }
        res = nls(~ProfileSSE.DDE(pars, allpars, times, data, coefs,
            lik, proc, in.meth, control.in, active, basisvals = basisvals, fdobj0=fdobj0),
        data = list(allpars = pars, basisvals = basisvals, fdobj0 = fdobj0,
            times = times, data = data, coefs = ncoefs, lik = lik,
            proc = proc, in.meth = in.meth, control.in = control.in,
            active = active), start = list(pars = pars[active]),
            trace = control.out$trace, control = control.out)
        apars = res$m$getPars()
        g = res$m$gradient()
        resid = res$m$resid()
        if (file.exists("curcoefs.tmp"))
            ncoefs = as.matrix(read.table(file = "curcoefs.tmp"))
        else ncoefs = coefs
    }
    names(apars) = aparamnames
    pars[active] = apars
    ncoefs = as.matrix(ncoefs)
    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
    if (file.exists("curcoefs.tmp")) {
        file.remove("curcoefs.tmp")
    }
    if (file.exists("optcoefs.tmp")) {
        file.remove("optcoefs.tmp")
    }
    if (file.exists("counter.tmp")) {
        file.remove("counter.tmp")
    }
    if (!is.null(fd.obj)) {
        ncoefs = array(ncoefs, c(nrow(ncoefs)/dims[2], dims[2],
            dims[3]))
        fd.obj = fd(ncoefs, fd.obj$basis)
        return(list(pars = pars, fd = fd.obj, lik = lik, proc = proc,
            outer.result = res))
    }
    else {
        return(list(pars = pars, coefs = ncoefs, lik = lik, proc = proc,
            outer.result = res, data = data, times = times))
    }
}

checkweights  <- function(weights, whichrows, diffs)
{
    if(is.vector(diffs)){
        diffs <- matrix(diffs, length(diffs), 1)
    }
    if (is.null(whichrows)) {
        whichrows = 1:nrow(diffs)
    }
    if (is.null(weights)) {
        return(matrix(1, nrow(diffs), ncol(diffs)))
    }
    else if (prod(dim(weights[whichrows, ]) == dim(diffs))) {

        if(dim(weights)[2] == 1){
            return(matrix(weights[whichrows, ],ncol = 1))
        }

        return(weights[whichrows, ])
    }
    else if (length(weights) == ncol(diffs)) {
        return(matrix(weights, nrow(diffs), ncol(diffs), byrow = TRUE))
    }
    else if (length(weights[whichrows]) == nrow(diffs)) {
        return(matrix(weights, nrow(diffs), ncol(diffs), byrow = FALSE))
    }
    else if (ncol(weights[whichrows, ]) > ncol(diffs) & nrow(weights) >
        nrow(diffs)) {
        warning("Dimension of weights does not match that of data")
        return(weights[whichrows[1:nrow(diffs)], 1:ncol(diffs)])
    }
    else {
        stop("Dimension of weights does not match that of data")
    }
}


Smooth.LS.DDE <- function(fn, data, times, pars, coefs = NULL, coefs0 = NULL,  basisvals = NULL, basisvals0 = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", control.in = list(), eps = 1e-06, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, tauMax = NULL, delay = NULL)
{
    if(is.null(tauMax)) tauMax <- min(times) + 1/3 * (range(times)[2]- range(times)[1])
    dims = dim(data)
    ## if(dims[2] == 1){
    ##     data <- matrix(data[times > pars["tau"],], ncol = 1)
    ## }else
    ##     data <- data[times > pars["tau"],]
    ## coefs <- coefs[times > pars["tau"]]
    ## times <- times[times > pars["tau"]]
    ## Prepare to calculate the delay:
    ## norder <- basisvals$nbasis - length(basisvals$params)
    ## nbasis.d <- length(knots.d) + norder - 2
    ## range.d <- range(knots.d)

    ## Create fd objects to fit y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"

    ## basisvalsTau <- basisvals
    ## basisvalsTau$rangeval[1] <- pars["tau"]
    ## basisvalsTau$params <- basisvals$params[ basisvals$params > pars["tau"]]
    ## basisvalsTau$nbasis <-  basisvals$nbasis - length(basisvals$params) +
    ##     length(basisvalsTau$params)
    ## basisvalsTau$names <- basisvals$names[1:basisvalsTau$nbasis]
    profile.obj <- LS.setup(pars, coefs, fn, basisvals=basisvals, lambda,
        fd.obj, more, data, weights, times, quadrature, eps = 1e-06,
        posproc, poslik, discrete, names = names, sparse = sparse,
        likfn = make.id(), likmore = NULL)
    lik = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times

    ##################################################
    ## Added delay data and functions
    ##################################################
    delayLikObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times, pars["tau"])
    delayProcObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, pars["tau"])
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$dy.d <- delayLikObj$dy.d
    proc$more$more$dy.d <- delayProcObj$dy.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    lik$more$more$dbvals.d <- delayLikObj$dbvals.d
    proc$more$more$dbvals.d <- delayProcObj$dbvals.d
    proc$dfdc <- delay$dfdc
    proc$d2fdc2 <- delay$d2fdc2.DDE
    proc$d2fdcdp <- delay$d2fdcdp.DDE
    proc$more$delay <- delay
    ## No need to use:
    ## proc$more$dfdtau <- dfdtau.DDE
    ## proc$more$d2fxdtau <- d2fxdtau.DDE
    ## proc$more$d2fdx.ddtau <- d2fdx.ddtau.DDE

    ## fdnames <- list(NULL, NULL, NULL)
    ## fdnames[[2]] <- attr(coefs0, "dimnames")[[2]]
    ## fdobj <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    ## attr(fdobj, "class") <- "fd"
    ## lik$more$more$y.d <- eval.fd(times - pars["tau"], fdobj)
    ## proc$more$more$y.d <- eval.fd(proc$more$qpts - pars["tau"], fdobj)

    Ires = inneropt.DDE(data, times, pars, coefs, lik, proc, in.meth,
        control.in, basisvals=basisvals, fdobj0 = fdobj0)
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

SplineCoefsErr.DDE <- function(coefs, times, data, lik, proc, pars, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = pars["tau"])
    delayLikObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = times, tau = pars["tau"])
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$dy.d <- delayLikObj$dy.d
    proc$more$more$dy.d <- delayProcObj$dy.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    lik$more$more$dbvals.d <- delayLikObj$dbvals.d
    proc$more$more$dbvals.d <- delayProcObj$dbvals.d

    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2,
        proc$bvals, pars, proc$more)
    if (!is.null(proc$report)) {
        print(f)
    }
    return(sgn * f)
}

SplineCoefsDC.DDE <- function(coefs, times, data, lik, proc, pars, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = pars["tau"])
    delayLikObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = times, tau = pars["tau"])
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$dy.d <- delayLikObj$dy.d
    proc$more$more$dy.d <- delayProcObj$dy.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    lik$more$more$dbvals.d <- delayLikObj$dbvals.d
    proc$more$more$dbvals.d <- delayProcObj$dbvals.d
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals,
        pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars,
        proc$more)
    g = as.vector(g)
    return(sgn * g)
}

SplineCoefsDC2.DDE <- function(coefs, times, data, lik, proc, pars, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = pars["tau"])
    delayLikObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = times, tau = pars["tau"])
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$dy.d <- delayLikObj$dy.d
    proc$more$more$dy.d <- delayProcObj$dy.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    lik$more$more$dbvals.d <- delayLikObj$dbvals.d
    proc$more$more$dbvals.d <- delayProcObj$dbvals.d
    result = as.matrix(SplineCoefsDC2sparse(coefs, times, data,
    lik, proc, pars, sgn))
    return(result)
}

inneropt.DDE <- function(data, times, pars, coefs, lik, proc,
                         in.meth = "nlminb", control.in = list(),
                         basisvals, fdobj0)
{
    check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    if (in.meth == "optim") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$maxit)) {
            control.in$maxit = 1000
        }
        if (is.null(control.in$reltol)){
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
        res = optim(coefs, SplineCoefsErr.DDE, gr = SplineCoefsDC.DDE,
        hessian = control.in$reportHessian, control = control.in,
        times = times, data = data, lik = lik, proc = proc,
        pars = pars, method = imeth, basisvals = basisvals, fdobj0 = fdobj0)
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
            Hessian = SplineCoefsDC2.DDE
        }
        else {
            Hessian = NULL
        }
        ## SplineCoefsErr do not need to be changed.
        ##
        res = nlminb(coefs, SplineCoefsErr.DDE, gradient = SplineCoefsDC.DDE,
        hessian = Hessian,
        control = control.in, times = times,
            data = data, lik = lik, proc = proc, pars = pars,
            basisvals = basisvals, fdobj0 = fdobj0)
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

## Added inputs:
## tauIndex, basisObj, fdobj0
ProfileSSE.AllPar.DDE <- function(pars, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1,
                             basisvals, fdobj0, tauMax = NULL)
{
    ## Squared Error outer criterion

    ##    coefs = as.vector(coefs)  # First run the inner optimization
    ##
    f1 = SplineCoefsErr.DDE(coefs,times,data,lik,proc,pars, basisvals = basisvals, fdobj0 = fdobj0)
    ## ?
    ## tau <- pars[proc$more$more$tauIndex]
    if(use.nls){
        # If using NLS, we need to keep track
        if(file.exists('curcoefs.tmp')){
          altcoefs = as.matrix(read.table('curcoefs.tmp'))
          if( !(length(altcoefs)==length(coefs)) ){
              stop(paste('Variables in curcoefs.tmp do not conform;',
                       'file exists from previous experiments?'))
          }
      }else {
          altcoefs = coefs
        }
        if(file.exists('counter.tmp')){
          counter = read.table('counter.tmp')
          niter = counter[nrow(counter),1]
        } else {
          counter = matrix(c(1,0,pars),1,length(pars)+2)
          niter = 0
      }

        f2 = SplineCoefsErr.DDE(altcoefs,times,data,lik,proc,pars, basisvals = basisvals, fdobj0 = fdobj0)

        if(f2 < f1){
            coefs = altcoefs
            f1 = f2
        }

        altdevals = as.matrix(lik$bvals %*% matrix(altcoefs,ncol(lik$bvals),
                              length(altcoefs)/ncol(lik$bvals)))
        colnames(altdevals) = proc$more$names
    }

    ## The coefficients for delayed times:
    if(is.null(tauMax)) tauMax <- min(times) + 1/3 * (range(times)[2]- range(times)[1])
    if(pars["tau"] > tauMax) pars["tau"] <- tauMax


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
    Ires = inneropt.DDE(data,times,pars,coefs,lik,proc, in.meth,control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ncoefs = Ires$coefs

    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj.d <- list(coefs = ncoefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    ## Can only deal with one "tau".
    delayLikObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times, pars["tau"])
    delayProcObj <- delay.fit(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, pars["tau"])
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$dy.d <- delayLikObj$dy.d
    proc$more$more$dy.d <- delayProcObj$dy.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    lik$more$more$dbvals.d <- delayLikObj$dbvals.d
    proc$more$more$dbvals.d <- delayProcObj$dbvals.d
    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names
    ## Squared errors: No need to change for DDE
    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
    f = as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more))*sqrt(weights))
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

    if(!is.null(lik$report)){ print(f) }
    f = sgn*f
    df = sgn*df

    ## Numerical derivative
    ## nParsDelay <- sum(proc$more$more$tauIndex)
    ## eps <- 1E-4
    ## for(i in (length(pars) - nParsDelay + 1):length(pars)){
    ##     pars.i1 <- pars.i2 <- pars
    ##     pars.i1[i] <- pars[i] + eps
    ##     pars.i2[i] <- pars[i] - eps
    ##     Ires.i1 <- inneropt.DDE(data,times,pars.i1,coefs,lik,proc, in.meth,control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ##     Ires.i2 <- inneropt.DDE(data,times,pars.i2,coefs,lik,proc, in.meth,control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ##     ncoefs.i1 <- Ires.i1$coefs
    ##     ncoefs.i2 <- Ires.i2$coefs
    ##     devals.i1 = as.matrix(lik$bvals%*%ncoefs.i1)
    ##     devals.i2 = as.matrix(lik$bvals%*%ncoefs.i2)
    ##     colnames(devals.i1) <- colnames(devals.i2) <- proc$more$names
    ##     f.i1 = as.vector( as.matrix(lik$more$fn(times, devals.i1, pars, lik$more$more))*sqrt(weights))
    ##     f.i2 = as.vector( as.matrix(lik$more$fn(times, devals.i2, pars, lik$more$more))*sqrt(weights))
    ##     df[,i] <- (f.i2 -  f.i1) / 2 / eps
    ## }

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




