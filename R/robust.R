##' Plot Diagnostics for an rlm Object
##'
##' @param object \code{\link[MASS]{rlm}} object
##' @param which a subset of the numbers \code{1:8}
##' @param ask logical. See \code{par("ask")}
##' @param label.pos  positiong of labels. See \code{plot.lm}
##' @param labels.id vector of labels. See \code{plot.lm}
##' @param cex.id magnification of point labels. See \code{plot.lm}
##' @param ... other parameters to be passed through to \code{plot.lm} function
##' @export
plot.rlm <- function(object, which=c(7L:8L), ask = prod(par("mfcol")) < length(which) && dev.interactive(), label.pos=c(4,2), labels.id = names(residuals(object)), cex.id=0.75, ...) {

    show <- rep(FALSE, 8)
    show[which] <- TRUE
    
    if(any(show[1L:6L])) {
        getS3method("plot", "lm")(object, which=which(show[1L:6L], label.pos=label.pos,
                                              labels.id=labels.id, cex.id=0.75, ...))
    }

    if (ask) {
	oask <- devAskNewPage(TRUE)
	on.exit(devAskNewPage(oask))
    }

    text.id <- function(x, y, ind, adj.x = TRUE) {
        labpos <-
            if(adj.x) label.pos[1+as.numeric(x > mean(range(x)))] else 3
        text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
             pos = labpos, offset = 0.25)
    }

    if (any(show[7L:8L])) {
        x <- model.matrix(object)[,-1]
        mcd <- MASS::cov.mcd(x)            
        mcd.dd <- sqrt(mahalanobis(x, mcd$center, mcd$cov))
    }
    if (show[7L]) {
        dev.hold()
        r <- with(object, residuals/s)
        plot(mcd.dd, r, xlab="Robust MCD Distance", ylab="Standardized Robust Residual")
        grid()
        abline(h=c(-3, 3), v=3)
        ids <- which(abs(r) > 3 | mcd.dd > 3)
        text.id(mcd.dd[ids], r[ids], ids)
        dev.flush()
    }
    if (show[8L]) {
        mahala.dd <- sqrt(mahalanobis(x, colMeans(x), cov(x)))

        dev.hold()
        plot(mahala.dd, mcd.dd, xlab="Mahalanobis Distance", ylab="Robust MCD Distance")
        abline(0, 1, lty=3)
        abline(h=3, lty=1)
        ids <- which(mcd.dd > 3)
        text.id(mahala.dd[ids], mcd.dd[ids], ids)
        dev.flush()
    }
    
    invisible()
}
