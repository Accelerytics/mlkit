% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/real.adaboost.R
\name{real.adaboost}
\alias{real.adaboost}
\title{Implementation of the Real AdaBoost Algorithm}
\usage{
real.adaboost(formula, data, classifier, n.iter)
}
\arguments{
\item{formula}{an object of class \code{\link[stats]{formula}}: a symbolic
description of the model to be fitted following the standard of
\code{\link[stats]{lm}}.}

\item{data}{an optional data frame, list or environment (or object coercible
by \code{\link[base]{as.data.frame}} to a data frame) containing the
variables in the model. If not found in data, the variables are taken from
environment (\code{formula}), typically the environment from which this
function is called.}

\item{classifier}{classifier to use in the algorithm.}

\item{n.iter}{number of iterations to apply.}
}
\value{
\code{real.adaboost} returns an object of \code{\link{class}}
\code{mlkit.class.fit}. An object of class \code{mlkit.class.fit} is a list
containing at least the following components:
\item{classifier}{a multiclass classifier function.}
\item{y.hat}{the in-sample classifications of the given data.}
}
\description{
Implementation of the Real AdaBoost algorithm based on
Friedman, J., Hastie, T., and Tibshirani, R. (2000)
}
