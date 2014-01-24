DSIRfn <- list()
DSIRfn$ode.fn <- function (t, y, parms)
{
    p = parms$p
    more = parms$more
    beta = more$beta.fun(t, p, more)
    tmpvec = beta * y[, "S"] * y[, "I"]
    r = y
    r[, "S"] = -tmpvec + p["mu"] - p["nu"] * y[, "S"]
    r[, "E"] = tmpvec - p["sigma"] * y[, "E"] - p["nu"] * y[,
        "E"]
    r[, "I"] = p["sigma"] * y[, "E"] - p["gamma"] * y[, "I"] -
        p["nu"] * y[, "I"]
    return(list(r))
}

DSIRfn$fn <- function (t, y, p, more)
{
    r = y
    yi.d <- more$yi.d
    r[, "S"] =  - p["beta"] * yi.d * y[, "S"]
    r[, "I"] = p["beta"] *  yi.d * y[, "S"] - p["gamma"] * y[, "I"]
    return(r)
}

DSIRfn$dfdx <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    yi.d <- more$yi.d

    r[, "S", "S"] = - p["beta"] * yi.d
    r[, "I", "S"] = p["beta"] * yi.d
    r[, "I", "I"] = -p["gamma"] * y[,"I"]
    return(r)
}

DSIRfn$dfdx.d <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(more$yi.d)))
    dimnames(r) = list(NULL, colnames(y), colnames(more$y.d))
    yi.d <- more$yi.d
    r[, "S", "I"] = - p["beta"] * y[,"S"]
    r[, "I", "I"] = p["beta"] * y[,"S"]
    return(r)
}



DSIRfn$dfdp <- function (t, y, p, more)
{
    yi.d <- more$yi.d
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    r[, "S", "beta"] = - yi.d * y[, "S"]
    r[, "I", "beta"] = yi.d * y[, "S"]
    r[, "I", "gamma"] = - y[, "I"]
    return(r)
}

DSIRfn$d2fdx2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

DSIRfn$d2fdxdp <- function (t, y, p, more)
{
    beta = more$beta.fun(t, p, more)
    dbetadp = more$beta.dfdp(t, p, more)
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    tmpdiag1 = diag(p["i"] + y[, "I"]) %*% dbetadp
    tmpdiag2 = diag(y[, "S"]) %*% dbetadp
    r[, "S", "S", more$beta.ind] = -tmpdiag1
    r[, "S", "I", more$beta.ind] = -tmpdiag2
    r[, "E", "S", more$beta.ind] = tmpdiag1
    r[, "E", "I", more$beta.ind] = tmpdiag2
    r[, "S", "S", "i"] = -beta
    r[, "E", "S", "i"] = beta
    r[, "S", "S", "nu"] = -1
    r[, "E", "E", "nu"] = -1
    r[, "I", "I", "nu"] = -1
    r[, "E", "E", "sigma"] = -1
    r[, "I", "E", "sigma"] = 1
    r[, "I", "I", "gamma"] = -1
    return(r)
}
