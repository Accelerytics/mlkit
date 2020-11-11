#' Fitting Linear Models using MM algorithm
#'
#' @description Implementation of the MM algorithm solver for a linear
#' regression model.
#'
#' @param formula an object of class \code{\link[stats]{formula}}: a symbolic
#' description of the model to be fitted following the standard of
#' \code{\link[stats]{lm}}.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link[base]{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in data, the variables are taken from
#' environment (\code{formula}), typically the environment from which this
#' function is called.
#' @param intercept optional boolean indicating whether to fit an intercept. If
#' \code{TRUE}, \code{standardize} is ignored. Default is \code{FALSE}.
#' @param standardize optional boolean indicating whether to return results for
#' standardized data. If \code{intercept} is \code{TRUE}, this argument is
#' ignored. Default is \code{FALSE}.
#' \code{TRUE}, \code{standardize} is ignored. Default is \code{FALSE}.
#' @param beta.init optional initial beta parameters to use in the MM
#' algorithm. Default is \code{NULL}.
#' @param beta.tol optional absolute tolerance for rounding down parameter
#' standardized estimates. If the absolute value of a parameter estimate in the
#' standardized model is smaller than \code{beta.tol}, it is rounded down to
#' zero. Default is \code{0}, that is, no rounding.
#' @param loss.tol optional convergence tolerance on the elastic net loss in
#' the MM algorithm. Default is \code{1e-6}.
#' @param seed optional seed. Default is \code{NULL}.
#' @param verbose optional number indicating per how many iterations the
#' estimation progress is displayed. Default is \code{0}, that is, no progress
#' updates.
#'
#' @return \code{mm.lm} returns an object of \code{\link{class}}
#' \code{mlfit}. An object of class \code{mlfit} is a list containing at
#' least the following components:
#' \item{coefficients}{a named vector of optimal coefficients.}
#' \item{loss}{residual sum of squares for optimal
#' coefficients.}
#' \item{r2}{coefficient of determination for optimal coefficients.}
#' \item{adj.r2}{adjusted coefficient of determination for optimal
#' coefficients.}
#' @export
#'
mm.lm = function(formula, data, intercept=F, standardize=F, beta.init=NULL,
  beta.tol=0, loss.tol=1e-6, seed=NULL, verbose=0) {

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

  # Define constants
  Xt.X = crossprod(x); inv.lambda = 1 / eigen(Xt.X)$values[1]
  inv.lambda.Xt.y = inv.lambda * crossprod(x, y); set.seed(seed)

  # Define helper functions for computing loss and displaying results
  loss = function(beta) rss(beta, x, y)
  pline = function(i, o, n, d)
    list(c('Iteration', i), c('Loss.old', o), c('Loss.new', n), c('Delta', d))

  # Choose some inital beta_0 and compute initial loss
  b.new = initialize.beta(beta.init, x); l.new = loss(b.new)

  # Update iteration and replace old parameters by previous until convergence
  i = 0L; while (TRUE) { i = i + 1L; l.old = l.new

    # Update parameters, loss and delta
    b.new = b.new - inv.lambda * Xt.X %*% b.new + inv.lambda.Xt.y
    l.new = rss(b.new); diff = l.old - l.new

    # Display progress if verbose
    if (verbose & (i %% verbose == 0))
      progress.str(pline(i, l.new, l.old, diff))

    # Break if improvement smaller than tol, that is, sufficient convergence
    if (diff / l.old < loss.tol) break
  }

  # Ensure information of last iteration is displayed
  if (verbose & (i %% verbose)) progress.str(pline(i, l.new, l.old, diff))

  # Force elements smaller than beta.tol to zero
  b.new[abs(b.new) < beta.tol] = 0

  # Descale beta if necessary
  if (intercept | standardize) beta = c(0, b.new)
  else beta = descale.beta(b.new, x.mean, x.sd, y.mean, y.sd)

  # Return mlfit object
  res = list('coefficients'=beta, 'alpha'=0, 'lambda'=0, 'loss'=loss(b.new),
    'R2'=r2(b.new, x, y), 'adj.R2'=adj.r2(b.new, x, y))
  class(res) = 'mlfit'
  return(res)
}
