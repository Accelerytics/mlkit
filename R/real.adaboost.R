#' Implementation of the Real AdaBoost Algorithm
#'
#' @description Implementation of the Real AdaBoost algorithm based on
#' Friedman, J., Hastie, T., and Tibshirani, R. (2000)
#'
#' @param formula an object of class \code{\link[stats]{formula}}: a symbolic
#' description of the model to be fitted following the standard of
#' \code{\link[stats]{lm}}.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link[base]{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in data, the variables are taken from
#' environment (\code{formula}), typically the environment from which this
#' function is called.
#' @param classifier classifier to use in the algorithm.
#' @param n.iter number of iterations to apply.
#'
#' @return \code{real.adaboost} returns an object of \code{\link{class}}
#' \code{mlkit.class.fit}. An object of class \code{mlkit.class.fit} is a list
#' containing at least the following components:
#' \item{classifier}{a multiclass classifier function.}
#' \item{y.hat}{the in-sample classifications of the given data.}
#'
#' @export
#'
real.adaboost = function(formula, data, classifier, n.iter) {

  # Extract dependent variable and explanatory variables
  y = as.vector((data[, all.vars(formula)[1]]))
  x = cbind(1, stats::model.matrix(formula, data))

  # Ensure y is correctly formatted and non-singular
  if (setequal(unique(y), c(0L, 1L))) y[y == 0L] = -1L
  if (!setequal(unique(y), c(-1L, 1L)))
    stop("Remap outcomes to +1 and -1 and ensure both are present.")

  # Step 1: Initialize variables
  w = rep(1 / n.obs, nrow=n.obs); f = list(); p.func = list()

  # Step 2: Apply n.iter updates
  for (m in 1:n.iter) {

    # a: Fit classifier to obtain class probabilities
    p.func[[m]] = function(x.new)
      predict(classifier(x=x, y=y, weights=w), x.new)

    # b: Specify f[[m]]
    f[[m]] = function(x) {
      p = p.func[[m]](x)
      return(0.5 * log(p / (1 - p)))
    }

    # c: Update and normalize weights
    w = w * exp(-y * f[[m]](x)); w = w / sum(w)
  }

  # Step 3: Construct classifier function to return
  classifier = function(x) sign(sum(sapply(1:n.iter, function(m) f[[m]](x))))

  # Return mlfit.class.fit object
  res = list(classifier=classifier, y.hat=classifier(x))
  class(res) = 'mlkit.class.fit'
  return(res)
}
