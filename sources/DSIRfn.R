DSIRfn <- list()

DSIRfn$fn <- function (t, y, p, more)
{
    r = y
    yi.d <- more$y.d[,2]
    r[, "S"] =  - p["beta"] * yi.d * y[, "S"] + (1 - y[,"S"]) * (0.25 * (sin(t/2) + 1) + 0.1)
    r[, "I"] = p["beta"] *  yi.d * y[, "S"] - (p["gamma"] + (0.25 * (sin(t/2) + 1) + 0.1)) * y[, "I"]
    return(r)
}

DSIRfn$dfdx <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    yi.d <- more$y.d[,2]
    r[, "S", "S"] = - p["beta"] * yi.d
    r[, "I", "S"] = p["beta"] * yi.d
    r[, "I", "I"] = -p["gamma"] - (0.25 * (sin(t/2) + 1) + 0.1)
    return(r)
}

DSIRfn$dfdx.d <- function (t, y, p, more)
{
    yi.d <- more$y.d[,2]
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[, "S", "I"] = - p["beta"] * y[,"S"]
    r[, "I", "I"] = p["beta"] * y[,"S"]
    return(r)
}



DSIRfn$dfdp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,2]
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
    yi.d <- more$y.d[,2]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "S", "S", "beta"] = - yi.d
    r[, "I", "S", "beta"] = yi.d
    r[, "I", "I", "gamma"] = - 1
    return(r)
}


DSIRfn$d2fdx.ddp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,2]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "S", "I", "beta"] = - y[,"S"]
    r[, "I", "I", "beta"] =  y[,"S"]
    return(r)
}

DSIRfn$d2fdxdx.d <- function (t, y, p, more)
{
    yi.d <- more$y.d[,2]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    r[, "S", "S", "I"] = - p["beta"]
    r[, "I", "S", "I"] = p["beta"]
    return(r)
}


DSIRfn$d2fdx.d2 <- function (t, y, p, more)
{
    yi.d <- more$y.d[,2]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

