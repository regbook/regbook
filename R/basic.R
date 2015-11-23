#' Fancy Histogram
#'
#' @param x a numeric vector
#' @param breaks see \link[graphics]{hist}
#' @examples
#' with(subset(hweight, gender=="M"), histf(height))
#' with(subset(hweight, gender=="M"), {
#'    hist(height, probability=TRUE, breaks=sqrt(length(height)), col="lightblue")
#'    lines(density(height))
#'    xs <- seq(min(height), max(height), length.out=100)
#'    lines(xs, dnorm(xs))
#' })
#' @export
histf <- function(x, breaks=sqrt(length(x)), length.out=100, rug=FALSE, boxplot=FALSE, fill="lightblue",
                  main=paste("Histogram of", xname),
                  densityFunction=function(y) dnorm(y, mean(x), sd(x)),
                  densityName="normal", xrange=extendrange(c(min(x), max(x))), 
                  ...) {
    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    xname <- deparse(substitute(x))

    hist.and.density <- function(...) {
        hist(x, probability=TRUE, breaks=breaks, xlim=xrange, col=fill, main=main, ...)
        lines(density(x), col=1, lty=3, ...)
        xgrid <- seq(xrange[1], xrange[2], length.out=length.out)
        lines(xgrid, densityFunction(xgrid), lty=1, col=2)
        legend("topright", legend=c("kernel", densityName), lty=c(3, 1), col=1:2)
        if(rug) rug(x)
    }
    if(boxplot) {
        layout(matrix(c(1, 2), 2, 1), heights = c(3, 1))
        par(mar=c(0, 4, 4, 2))    
        hist.and.density(xaxt="n", xlab="")
        par(mar=c(5, 4, 0, 2))
        boxplot(x, horizontal=TRUE, frame=FALSE, ylim=xrange, xlab=xname)
    } else {
        hist.and.density(...)
    }
}

#'The Bivariate Normal Distribution
#'
#' @param x a 2-dimensional vector or a 2 by n matrix
#' @param mean a mean vector
#' @param sd a standard deviation vector
#' @param rho correlation coefficient
#' @examples
#' x <- seq(-3,3,length.out=50)
#' y <- seq(-3,3,length.out=50)
#' z <- outer(x, y, function(x, y) dbvnorm(cbind(x, y), rho=0.5))
#' persp(x,y,z,theta=5,phi=20,col="lightblue", expand=0.5, zlab="density")
#' contour(x,y,z, xlab="x", ylab="y")
#' @export
dbvnorm <- function(x, mean=c(0, 0), sd=c(1, 1), rho=0) {
    if (is.vector(x)) x <- matrix(x, ncol=2)
    z <- (t(x) - mean)/sd
    s <- 1 - rho^2
    k <- 1/((2*pi)*prod(sd)*sqrt(s))
    A <- z[1,]^2 + z[2,]^2 - 2*rho*z[1,]*z[2,]
    k * exp(-A/(2*s))
}



