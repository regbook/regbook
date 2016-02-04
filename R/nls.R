#' Wald confidence interval
#' @export
waldint <- function(obj, ...) {
    UseMethod("waldint")
}


#' @rdname waldint
#' @export
waldint.nls <- function(obj, parm, level = 0.95, ...) {
    cf <- coef(obj)
    pnames <- names(cf)
    if(missing(parm)) parm <- seq_along(pnames)
    if(is.numeric(parm)) parm <- pnames[parm]
    
    b <- cf[parm]
    s <- sqrt(diag(vcov(obj)))[parm]

    a <- (1-level)/2

    lwb <- b + s*qt(a, df.residual(obj))
    upb <- b + s*qt(1-a, df.residual(obj))
    ans <- cbind(lwb, upb)
    colnames(ans) <- paste(format(c(a, 1-a)*100, 3), "%", sep="")
    ans
}
