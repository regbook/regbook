% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/robust.R
\name{plot.rlm}
\alias{plot.rlm}
\title{Plot Diagnostics for an rlm Object}
\usage{
\method{plot}{rlm}(object, which = c(7L:8L), ask = prod(par("mfcol")) <
  length(which) && dev.interactive(), label.pos = c(4, 2),
  labels.id = names(residuals(object)), cex.id = 0.75, ...)
}
\arguments{
\item{object}{\code{\link[MASS]{rlm}} object}

\item{which}{a subset of the numbers \code{1:8}}

\item{ask}{logical. See \code{par("ask")}}

\item{label.pos}{positiong of labels. See \code{plot.lm}}

\item{labels.id}{vector of labels. See \code{plot.lm}}

\item{cex.id}{magnification of point labels. See \code{plot.lm}}

\item{...}{other parameters to be passed through to \code{plot.lm} function}
}
\description{
Plot Diagnostics for an rlm Object
}

