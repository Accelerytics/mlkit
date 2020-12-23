#' Implementation of the Multiclass LogitBoost Algorithm
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
#'
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
  dep.var = all.vars(formula)[1]
  y = as.factor(df$train[[dep.var]])
  y.star = stats::model.matrix(~ 0 + ., as.data.frame(y))
  x = cbind(1, stats::model.matrix(formula, data))

  # Ensure y is correctly formatted and non-singular
  n.obs = nrow(y.star); n.class = ncol(y.star)
  if (n.class == 1 | n.class == n.obs)
    stop("Ensure dependent variable contains between 1 and n.obs categories.")

  # Step 1: Initialize variables
  w = matrix(1 / n.obs, nrow=n.obs, ncol=n.class); f = list()
  z = matrix(nrow=n.obs, ncol=n.class); upt.w = (n.class - 1) / n.class
  p = matrix(1 / n.class, nrow=n.obs, ncol=n.class)

  # Step 2: Apply n.iter updates to find MLE
  for (m in 1:n.iter) {

    # a.i: Compute working responses and weights in all classes
    w = p * (1 - p); z = (y.star - p) / w

    # a.ii: Fit WLS regression
    f[[m]] = lapply(1:n.class,
      function(j) glm.fit(x=x, y=z[, j], weights=w[, j])$coefficients)

    # b: Boosting update
    f[[m]] = lapply(1:n.class,
      function(j) upt.w * (f[[m]][[j]] - Reduce('+', f[[m]]) / n.class))
    f.func = function(x) sapply(
      1:n.class,
      function(j) x %*% Reduce('+', lapply(1:m, function(k) f[[k]][[j]]))
    )

    # c: Update probabilities
    tmp = exp(f.func(x)); p = tmp / sum(tmp)
  }

  # Step 3: Construct classifier function to return
  classes = unique(y)
  classifier = function(x) classes[apply(f.func(x), 1, which.max)]

  # Return mlfit.class.fit object
  res = list(classifier=classifier, y.hat=classifier(x))
  class(res) = 'mlkit.class.fit'
  return(res)
}
