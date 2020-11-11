#' Fitting Linear Models with Ridge Penalty
#'
#' @description Implementation of the analytical solution for a linear
#' regression model with a Ridge penalty term.
#'
#' @param formula an object of class \code{\link[stats]{formula}}: a symbolic
#' description of the model to be fitted following the standard of
#' \code{\link[stats]{lm}}.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link[base]{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in data, the variables are taken from
#' environment (\code{formula}), typically the environment from which this
#' function is called.
#' @param lambda penalty term scaling hyperparameter.
#' @param intercept optional boolean indicating whether to fit an intercept. If
#' \code{TRUE}, \code{standardize} is ignored. Default is \code{FALSE}.
#' @param standardize optional boolean indicating whether to return results for
#' standardized data. If \code{intercept} is \code{TRUE}, this argument is
#' ignored. Default is \code{FALSE}.
#' \code{TRUE}, \code{standardize} is ignored. Default is \code{FALSE}.
#' @param beta.tol optional absolute tolerance for rounding down parameter
#' standardized estimates. If the absolute value of a parameter estimate in the
#' standardized model is smaller than \code{beta.tol}, it is rounded down to
#' zero. Default is \code{0}, that is, no rounding.
#'
#' @return \code{ridge.lm} returns an object of \code{\link{class}}
#' \code{mlfit}. An object of class \code{mlfit} is a list containing at
#' least the following components:
#' \item{coefficients}{a named vector of optimal coefficients.}
#' \item{loss}{residual sum of squares plus ridge loss for optimal
#' coefficients.}
#' \item{r2}{coefficient of determination for optimal coefficients.}
#' \item{adj.r2}{adjusted coefficient of determination for optimal
#' coefficients.}
#' @export
#'
ridge.lm = function(formula, data, lambda, intercept=F, standardize=F,
  beta.tol=0) {

  # Extract dependent variable and explanatory variables
  y = data.matrix(data[, all.vars(formula)[1]]);
  x = stats::model.matrix(formula, data)

  # Store scaling parameters
  if (!standardize) {
    x.mean = colMeans(x); x.sd = apply(x, 2, stats::sd)
    y.mean = mean(y); y.sd = stats::sd(y)
  }

  # Add intercept or standardize data if necessary
  x = create.x(x, intercept, standardize)
  y = create.y(y, intercept, standardize)

  # Derive beta estimate
  b.new = solve(crossprod(x) + nrow(x) * lambda / 2 * diag(ncol(x)),
    crossprod(x, y))

  # Set elements smaller than beta.tol to zero
  b.new[abs(b.new) < beta.tol] = 0

  # Descale beta if necessary
  if (intercept | !standardize) beta = c(0, b.new)
  else beta = descale.beta(b.new, x.mean, x.sd, y.mean, y.sd)

  # Return mlfit object
  res = list('coefficients'=beta, 'alpha'=0, 'lambda'=lambda,
    'loss'=elastic.net.loss(b.new, x, y, intercept, 0, lambda),
    'R2'=r2(b.new, x, y), 'adj.R2'=adj.r2(b.new, x, y))
  class(res) = 'mlfit'
  return(res)
}
