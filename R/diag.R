#' Variance Inflation Factor and Collinearity Diagnostics
#' @param mod an R object
#' @param ... further arguments pass to or from other methods
#' @export
vif <- function(mod, ...) {
    UseMethod("vif")
}

#' @rdname vif
#' @export
vif.lm <- function(mod, ...) {
    mat <- model.matrix(mod)
    S <- vcov(mod)
    if(colnames(S)[1] == "(Intercept)") {
        mat <- mat[,-1]
        S <- S[-1, -1]
    }
    R <- cov2cor(S)
    k <- ncol(R)
    VIF <- sapply(1:k, function(i) det(R[-i, -i])) / det(R)
    names(VIF) <- colnames(R)

    ans <- list(VIF=VIF, design.matrix=mat)
    class(ans) <- "vif"
    ans
}

#' @rdname vif
#' @export
print.vif <- function(obj, ...) {
    print(obj$VIF)
}

#' @rdname vif
#' @export
summary.vif <- function(obj, ...) {

    R <- cor(obj$design.matrix)
    varnames <- colnames(R)
    eig <- eigen(R)
    eigval <- eig$values
    names(eigval) <- varnames
    
    condind <- sqrt(eig$values[1] / eig$values)

    propvar <- prop.table(t(eig$vectors^2)/eig$values, margin=2)
    colnames(propvar) <- varnames

    ans <- list(VIF=obj$VIF, correlation=R, eigenvalues=eig$values, condition.index=condind, proportion.of.variation=propvar)
    class(ans) <- "summary.vif"
    ans
}

#' @rdname vif
#' @export 
print.summary.vif <- function(obj, digits=max(3L, getOption("digits") - 3L), ...) {
    cat("\nVIF:\n")
    print.default(obj$VIF, digits = digits, ...)
    cat("\nProportion of Variation:\n")
    mat <- cbind(Eigenvalues=obj$eigenvalues, "Condition Index"=obj$condition.index, obj$proportion.of.variation)
    rownames(mat) <- 1:nrow(mat)
    print.default(mat)
    
}



##' Extract Standardized Coefficients
##'
##' @param mod 
##' @param ... 
##' @export
stdcoef <- function(mod, ...) {
    UseMethod("stdcoef")
}

#' @rdname stdcoef
##' @export
stdcoef.lm <- function(mod, ...) {
  X <- model.matrix(mod)
  y <- model.response(mod$model)
  Sjj <- diag(t(X) %*% X) - nrow(X) * colMeans(X)^2
  SYY <- sum(y^2) - length(y) * mean(y)^2
  sqrt(Sjj/SYY) * coef(mod)
}

##' Add Standardized Coefficients
##'
##' @param mod 
##' @param ...
##' @export
lmbeta <- function(mod, ...) {
    UseMethod("lmbeta")
}

##' @rdname lmbeta
##' @export
lmbeta.lm <- function(mod, ...) {
  mod$standardized.coefficients <- stdcoef(mod)
  class(mod) <- c("lmbeta", class(mod))
  mod
}

##' @rdname lmbeta
##' @export
print.lmbeta <- function(obj, ...) {
  cat("\nCall:\n")
  print(obj$call, ...)
  cat("\nStandardized Coefficients:\n")
  print(obj$standardized.coefficients, ...)
  cat("\nCoefficients:\n")
  print(obj$coefficients)
  cat("\n")
}

##' @rdname lmbeta
##' @export
summary.lmbeta <- function(obj, ...) {
  ans <- summary.lm(obj)
  ans$coefficients <- cbind("Std. Coef"=obj$standardized.coefficients, ans$coefficients)
  ans
}
