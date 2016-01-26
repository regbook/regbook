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
    cat("\nVariance Proportion:\n")
    mat <- cbind(Eigenvalues=obj$eigenvalues, "Cond.Index"=obj$condition.index, obj$proportion.of.variation)
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
  ans$coefficients <- cbind("BETA"=obj$standardized.coefficients, ans$coefficients)
  ans
}

##' Mean and Variance of Data Subsets
##'
##' Splits the data into subsets according to the X levels, compute mean and
##' variance of the response for each.
##' 
##' @param formula a formula, such as \code{y ~ x}, where the \code{y} variable are numeric data to be split into groups according to the grouping \code{x} variables.
##' @param data a data frame
##' @examples
##' meanvar(Y ~ X, restaurant)
##' plot(meanvar(Y ~ X, restaurant), sd ~ mean + I(mean^2))
##' @export
##' 
meanvar <- function(formula, data, ..., na.action = na.omit) {
  agg <- aggregate(formula, data=data,
                   FUN=function(x) c(mean=mean(x), var=var(x)),
                   ..., na.action = na.action)
  nc <- ncol(agg)
  xlevels <- agg[1:(nc-1)]
  m <- agg[[nc]][,1]
  v <- agg[[nc]][,2]

  ans <- list(xlevels=xlevels, mean=m, var=v, call=match.call())
  class(ans) <- c("meanvar")
  ans
}

##' @rdname meanvar
##' @export
as.data.frame.meanvar <- function(obj, ...) {
    data.frame(obj$xlevels, mean=obj$mean, var=obj$var)
}

##' @rdname meanvar
##' @export
print.meanvar <- function(obj, ...) {
  cat("\nCall:\n")
  print(obj$call)
  cat("\nMean and Variance:\n")
  print(data.frame(xlevels=obj$xlevels, mean=obj$mean, var=obj$var))
}

##' @rdname meanvar
##' @param obj an \code{meanvar} object
##' @param formula mean-vairance relationship, for example, \code{sd ~ mean}, \code{sd ~ mean + I(mean^2)}, \code{var ~ mean}, etc.
##' @export
plot.meanvar <- function(obj, formula,...) {
   if(missing(formula)) {
     plot(obj$mean, obj$var, xlab="Mean", ylab="Variance")
   } else {
     nw <- data.frame(mean=seq(min(obj$mean), max(obj$mean), length.out=100))
     if (formula[[2]] == "sd") {
       ylab <- "Standard Deviation"
       data <- data.frame(mean=obj$mean, sd=sqrt(obj$var))
     } else if(formula[[2]] == "var") {
       ylab <- "Variance"
       data <- data.frame(mean=obj$mean, var=obj$var)
     } else {
       stop("formula must be 'sd ~ ' or 'var ~ '")
     }
     plot(data[,1], data[,2], xlab="Mean", ylab=ylab)
     lines(nw$m, predict(lm(formula, data=data), newdata=nw))
  }
}

