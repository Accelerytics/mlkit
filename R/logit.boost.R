#' Implementation of the Binary LogitBoost Algorithm
#'
#' @description Implementation of the LogitBoost algorithm based on
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
#' @param n.iter number of iterations to apply.
#'
#' @return \code{logic.boost} returns an object of \code{\link{class}}
#' \code{mlkit.class.fit}. An object of class \code{mlkit.class.fit} is a list
#' containing at least the following components:
#' \item{classifier}{a multiclass classifier function.}
#' \item{y.hat}{the in-sample classifications of the given data.}
#'
#' @export
#'
logit.boost = function(formula, data, n.iter) {

  # Extract dependent variable and explanatory variables
  y = data[[all.vars(formula)[1]]]; n.obs = length(y)
  x = stats::model.matrix(formula, data)

  # Step 1: Initialize variables
  f = matrix(0, nrow=n.iter, ncol=ncol(x)); p = rep(0.5, n.obs)

  # Step 2: Apply n.iter updates to find MLE
  for (m in 1:n.iter) {

    # a.i: Compute working responses and weights in all classes
    w = pmax(p * (1 - p), 2 * .Machine$double.xmin)
    z = pmax(pmin(ifelse(y == 1, 1 / p, - 1 / (1 - p)), 4), -4)

    # a.ii: Fit WLS regression
    f[m, ] = glm.fit(x=x, y=z, weights=w)$coefficients

    # b: Boosting update
    f.func = function(x) 0.5 * x %*% colSums(f)

    # c: Update probabilities
    f.eval = f.func(x)
    p = exp(f.eval) / (exp(f.eval) + exp(-f.eval))
  }

  # Step 3: Construct classifier function to return
  classifier = function(x) sign(x %*% colSums(f))

  # Return mlfit.class.fit object
  res = list(classifier=classifier, y.hat=classifier(x))
  class(res) = 'mlkit.class.fit'
  return(res)
}
