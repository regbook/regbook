
#' Predicted Residuals
#'
#' Predicted Residuals
#' @param model an R object, typically by \code{lm}
#' @param infl influence structure
#' @param ... further arguments pass to or from other methods
#' @export
rpredict <- function(model, ...) UseMethod("rpredict")

#' @rdname rpredict
#' @export
rpredict.lm <- function(model, infl = lm.influence(model, do.coef=FALSE)) {
    res <- infl$wt.res / (1 - infl$hat)
    res[is.infinite(res)] <- NaN
    res
}

#' Predicted Residual Sum of Squares
#'
#' @export
press <- function(model, ...) UseMethod("press")

#' @rdname press
#' @export
press.lm <- function(model, ...) {
    pr <- rpredict(model)
    prss <- sum(pr^2)

    y <- model.response(model$model)
    sst <- sum((y - mean(y))^2)
    pr.r.sq <- 1 - prss / sst
    pr.r.sq <- ifelse(pr.r.sq < 0, 0, pr.r.sq)
    list(PRESS=prss, pred.r.squared=pr.r.sq)
}



##' Fancy Summary
##' 
##' @export
summaryf <- function(object, ...) UseMethod("summaryf")

##' Alternative summary for regsubsets object
##'
##' @param object a \code{\link[leaps]{regsubsets}} object
##' @param ...
##' @export
summaryf.regsubsets <- function(object, ...) {
    smry <- summary(object, ...)
    class(smry) <- "summaryf.regsubsets"
    smry
}

##' @rdname summaryf.regsubsets
##' @export
print.summaryf.regsubsets <- function(object, ...)  {
    with(object,
         print.data.frame(data.frame(outmat, rss, rsq, adjr2, cp, bic)))
}

##' Compute Predicted Residual Sum of Squares
##' 
##' @param object a \code{\link[leaps]{regsubsets}} object
##' @param ... 
##' @export
press.regsubsets <- function(object, ...) {
    vars <- all.vars(object$call)
    respname <- vars[1]

    dfname <- vars[length(vars)]
    x <- eval(as.name(dfname))

    inds <- summary(object)$which[,-1]
    cnames <- colnames(inds)

    n <- nrow(inds)
    PRESS <- numeric(n)
    pred.r.squared <- numeric(n)
    for (i in 1:n) {
        ind <- inds[i, ]
        fit <- lm(as.formula(paste(respname, paste(cnames[ind], collapse="+"), sep="~")), x)
        pr <- press(fit)
        PRESS[i] <- pr$PRESS
        pred.r.squared[i] <- pr$pred.r.squared
    }
    
    ans <- list(which=inds, PRESS=PRESS, pred.r.squared=pred.r.squared, outmat=summary(object)$outmat, obj=object)
    class(ans) <- "press.regsubsets"
    ans
}

##' @rdname press.regsubsets
##' @param ... 
##' @export
print.press.regsubsets <- function(object, ...) {
    with(object, print.data.frame(data.frame(outmat, PRESS, pred.r.squared)))
}
