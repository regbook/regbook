
#' Predicted Residuals
#'
#' Predicted Residuals
#' @param model an R object, typically by \code{lm}
#' @param infl influence structure
#' @param ... further arguments pass to or from other methods
#' @export
rpredict <- function(model, ...) UseMethod("rpredict")

#' @rdname rpredict
rpredict.lm <- function(model, infl = lm.influence(model, do.coef=FALSE)) {
    res <- infl$wt.res / (1 - infl$hat)
    res[is.infinite(res)] <- NaN
    res
}

#' Predicted Residual Sum of Squares
#' 
press <- function(object) {
    pr <- 
    prss <- sum(pr^2)

    y <- model.response(object$model)
    sst <- sum((y - mean(y))^2)
    pr.r.sq <- 1 - prss / sst
    list(PRESS=prss, pred.r.squared=pr.r.sq)
}