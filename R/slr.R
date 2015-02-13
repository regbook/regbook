#' Fitted Objectel Plot for Simple Linear Regression
#
#' @param object lm object. simple linear regression.
#' @param level tolerance/confidence level
#' @param legend legend location. See \code{\link[graphics]{legend}}
fitplot <- function(object, level=0.95,
                    legend=c("none", "topleft")) {
    mod <- object
    xname <- all.vars(delete.response(terms(mod)))
    if(length(xname) > 1) {
        stop("invalid model: only simple linear regression")
    }
    x <- model.matrix(mod)[, xname]
    
    ## compute confidence and prediction intervals
    xs <- with(aflength, seq(min(x), max(x), length.out=100))
    newdata <- data.frame(xs)
    colnames(newdata) <- xname
    ci <- predict(mod, newdata=newdata, interval="confidence", level=level)
    pr <- predict(mod, newdata=newdata, interval="prediction", level=level)

    plot(formula(mod), mod$model, ylim=range(c(ci, pr)))
    ##- confidence limits
    polygon(c(xs, rev(xs)), c(ci[,"lwr"], rev(ci[,"upr"])), col="#00000022", lty=0)
    ##- fitted line
    abline(mod, lwd=2, col="gray50")
    ##- prediction limits
    lines(xs, pr[,"lwr"], lty=2)
    lines(xs, pr[,"upr"], lty=2)
    
    ##- legend
    legend <- match.arg(legend)
    if(legend != "none") {
        legend(legend, legend=c("fit", paste(level*100, "% confidence limits", sep=""),
                              paste(level*100, "% prediction limits", sep="")),
               col=c("gray50", NA, 1), lwd=c(2, NA, 1), lty=c(1,NA,2), fill=c(NA, "#00000022", NA),
               border=c(NA, NA, NA), merge=TRUE)
    }
}
