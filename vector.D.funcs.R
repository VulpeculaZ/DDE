## y.d is added for delay
## vectorD.Fun <- vectorD$proc$more$fn <- function (times, y, p, more)
## {
##     y.d <- more$y.d
##     r <- y
##     r[, "S"] <- p["b"] * y.d[,"S"] * (1 - y[,"S"]) - p["a"] * y[,"S"]
##     return(r)
## }

## vectorD$proc$more$dfdx <- function(t, y, p, more){
##     y.d <- more$y.d
##     r = array(0, c(dim(y), dim(y)[2]))
##     dimnames(r) = list(NULL, colnames(y), colnames(y))
##     r[,"S", "S"] <- -p["b"] * y.d[,"S"] - p["a"]
##     return(r)
## }

## vectorD$proc$more$dfdp <- function(times, y, p, more)
## {
##     y.d <- more$y.d
##     r = array(0, c(dim(y), length(p)))
##     dimnames(r) = list(NULL, colnames(y), names(p))
##     r[, "V", "a"] <- - y[,"S"]
##     r[, "V", "b"] <- y.d[,"S"] * (1 - y[,"S"])
##     return(r)
## }

## ## Return an array of all 0's
## vectorD$proc$more$d2fdx2 <- function(times, y, p, more){
##     y.d <- more$y.d
##     r = array(0, c(dim(y), 1, 1))
##     dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
##     return(r)
## }

## vectorD$proc$more$d2fdxdp <- function (times, y, p, more)
## {
##     y.d <- more$y.d
##     r = array(0, c(dim(y), dim(y)[2], length(p)))
##     dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
##     r[,"S", "S", "a"] <- -1
##     r[,"S", "S", "b"] <- - y.d[,"S"]
## }

## ## New functions:
## vectorD$proc$more$dfdx.d <- function(times, y, p, more){
##     y.d <- more$y.d
##     r = array(0, c(dim(y), dim(y)[2]))
##     dimnames(r) = list(NULL, colnames(y), colnames(y.d))
##     r[,"S", "Delay"] <- p["b"] * (1 - y[,"S"])
##     return(r)
## }

## vectorD$proc$more$d2fdxdx.d <- function(times, x, p, more){
##     y.d <- more$y.d
##     r = array(-1, c(dim(y), 1, 1))
##     dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y.d))
##     return(r)
## }

## vectorD$proc$more$d2fdxdx.d <- function(times, x, p, more){
##     y.d <- more$y.d
##     r = array(0, c(dim(y), 1, 1))
##     dimnames(r) = list(NULL, colnames(y), colnames(y.d), colnames(y.d))
##     return(r)
## }

vectorFun <- list()

vectorFun$fn <- function (times, y, p, more)
{
    y.d <- more$y.d
    r <- y
    dimnames(r) <- dimnames(y) ## ??
    r[, "S"] <- (p["b"] + sin(times))* y.d[,"S"] * (1 - y[,"S"]) - p["a"] * y[,"S"]
    return(r)
}

##
vectorFun$fn.ode <- function(times, y, p)
{
    y.d <- more$y.d
    r <- y
    dimnames(r) <- dimnames(y) ## ??
    r[, "S"] <- (p["b"] + sin(times)) * y.d[,"S"] * (1 - y[,"S"]) - p["a"] * y[,"S"]
    return(r)
}

##
vectorFun$dfdx <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[,"S", "S"] <- -(p["b"] + sin(times)) * y.d[,"S"] - p["a"]
    return(r)
}

vectorFun$dfdp <- function(times, y, p, more)
{
    y.d <- more$y.d
    r = array(0, c(dim(y), 2))
    dimnames(r) = list(NULL, colnames(y), names(p)[c(1,2)])
    r[, "S", "a"] <- - y[,"S"]
    r[, "S", "b"] <- y.d[,"S"] * (1 - y[,"S"])
    return(r)
}


vectorFun$d2fdxdp <- function (times, y, p, more)
{
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2], 2))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p)[c(1,2)])
    r[,"S", "S", "a"] <- -1
    r[,"S", "S", "b"] <- - y.d[,"S"]
    return(r)
}

vectorFun$d2fdx.ddp <- function(times, y, p, more){
    y.d <- more$y.d
    r <- array(0, c(dim(y), dim(y)[2], 2))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p)[c(1:2)])
    r[,"S", "S", "a"] <- 0
    r[,"S", "S", "b"] <- 1 - y
    return(r)
}

vectorFun$d2fdx2 <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

vectorFun$dfdx.d <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y.d))
    r[,"S", "S"] <- (p["b"] + sin(times)) * (1 - y[,"S"])
    return(r)
}

vectorFun$d2fdx.ddx <- vectorFun$d2fdxdx.d <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(-(p["b"] + sin(times)), c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y.d))
    return(r)
}

vectorFun$d2fdx.d2 <- function(times, y, p, more){
    y.d <- more$y.d
    r = array(0, c(dim(y), 1, 1))
    dimnames(r) = list(NULL, colnames(y), colnames(y.d), colnames(y.d))
    return(r)
}

