##' Principal Component Regression
##'
##' @export
pcr <- function(formula, data, ncomp, cor=TRUE, ...) {
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)
    mt <- attr(mf, "terms")
    
    fit <- pcr.fit(x[,-1], y, ncomp=ncomp, cor=cor)
    class(fit$lm) <- "lm"
    fit$lm$terms <- mt 
    class(fit) <- "pcr"
    fit
}

##' Fitting Principal Component Regression Model
##' 
##' @param x explanatory matrix
##' @param y response vector
##' @param ncomp number of principal components
##' @param cor a logical value indicating whether the calculation should use
##' the correlation matrix or the covariance matrix. See \code{\link[stats]{princomp}}
##' @param ...
##' @export
pcr.fit <- function(x, y, ncomp, cor=TRUE, ...) {
    x <- as.matrix(x)
    y <- as.numeric(y)

    ## principal components
    pr <- princomp(x, cor=cor)
    pc <- pr$score

    ## fit principal component regression model
    fit <- lm.fit(cbind(1, pc[,1:ncomp]), y)

    list(ncomp=ncomp, princomp=pr, lm=fit)
}

##' @rdname pcr
##' @export
summary.pcr <- function(object, ...) {
    summary.lm(object$lm)
}

##' @rdname pcr
##' @export
print.pcr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    print(x$lm)
    ##print.default(format(x$coefficients, digits=digits), print.gap=2L, quote=FALSE)
}

##' Extract Model Coefficients
##' @export
coef.pcr <- function(x, type=c("original", "princomp")) {
    type <- match.arg(type)
    if(type == "original") {
        pr <- x$princomp
        fit <- x$lm
        ncomp <- x$ncomp
        b1 <- loadings(pr)[,1:ncomp] %*% coef(fit)[2:(ncomp+1)] / pr$scale
        names(b1) <- attr(fit$terms, "term.labels")
        b0 <- coef(fit)[1] - sum(b1 * pr$center)
        names(b0) <- "(Intercept)"
        c(b0, b1)
    } else if(type == "princomp") {
        coef(x$lm)
    }
}
