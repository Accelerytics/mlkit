#' Fitting Binary Support Vector Machine
#'
#' @description Implementation of the MM algorithm solver for the binary
#' support vector machine.
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
#' @param loss hinge loss to use in estimation. Default is \code{'qua'}, that
#' is, quadratic hinge loss.
#' @param huber.k hyperparameter for Huber hinge errors. Default is
#' \code{NULL}.
#' @param loss.tol optional convergence tolerance in the MM algorithm. Default
#' is \code{1e-6}.
#' @param v.init optional initial v parameters to use in the MM algorithm.
#' Default is \code{NULL}.
#' @param seed optional seed. Default is \code{NULL}.
#' @param verbose optional number indicating per how many iterations the
#' estimation progress is displayed. Default is \code{0}, that is, no progress
#' updates.
#'
#' @return \code{svm.bin} returns an object of \code{\link{class}}
#' \code{mlkit.bin.fit}. An object of class \code{mlkit.bin.fit} is a list
#' containing at least the following components:
#' \item{coefficients}{a named vector of optimal coefficients.}
#' \item{loss}{SVM loss including penalizing term on the parameters.}
#'
#' @export
#'
svm.bin = function(formula, data, lambda, loss='qua', huber.k=NULL,
  loss.tol=1e-6, v.init=NULL, seed=NULL, verbose=0) {

  # Hinge Huber loss with k equal to 1 equivalent to easier absolute hinge loss
  if (!is.null(huber.k))
    if ((loss == 'hub') & (huber.k == 1)) return(svm.bin(formula, data, lambda,
      'abs', loss.tol, v.init, seed, verbose))

  # Extract dependent variable and explanatory variables
  y = as.vector((data[, all.vars(formula)[1]]))
  x = cbind(1, stats::model.matrix(formula, data))

  # Ensure y is correctly formatted and non-singular
  if (setequal(unique(y), c(0L, 1L))) y[y == 0L] = -1L
  if (!setequal(unique(y), c(-1L, 1L)))
    stop("Remap outcomes to +1 and -1 and ensure both are present.")

  # Define constants
  n = nrow(x); POS.IND = (y == 1L); p = diag(c(0, rep(1, ncol(x) - 1)))

  # Define helper functions for computing specific expressions and loss
  h.loss = hinge.loss(x, loss, huber.k)
  upt.a = function(q) svm.upt.a(q, POS.IND, loss, huber.k)
  upt.b = function(a, q) svm.upt.b(a, q, POS.IND, loss, huber.k)

  # Updates of v depend on type of loss; use simplest for quadratic hinge loss
  if (loss == 'qua') {
    z = solve(t(x) %*% x + lambda * p) %*% t(x)
    update.v = function(a, x, lambda, b, v) return(z %*% b)
  } else {
    update.v = function(a, x, lambda, b, v) {
      tryCatch ({
        v.plus = solve(t(x) %*% a %*% x + lambda * p, t(x) %*% b)
        v + 2 * (v.plus - v)
      },
      error = function(e) {
        warning("Terminating early due to A being computationally singular.")
        v
      })
    }
  }
  t.loss = function(q, v) sum(unlist(sapply(q[!POS.IND & (q > -1)] + 1,
    h.loss))) + sum(unlist(sapply(1 - q[POS.IND & (q <= 1)], h.loss))) +
    lambda * sum(v[-1] ^ 2)
  pline = function(i, o, n, d)
    list(c('Iteration', i), c('Loss.old', o), c('Loss.new', n), c('Delta', d))

  # Choose some inital v_0 and compute initial loss
  set.seed(seed)
  v.new = initialize.beta(v.init, x); q = x %*% v.new; l.new = t.loss(q, v.new)

  # Update iteration and replace old parameters by previous until convergence
  i = 0L; while (TRUE) { i = i + 1L; l.old = l.new

    # Update parameters, loss and delta
    a = upt.a(q); v.new = update.v(a, x, lambda, upt.b(a, q), v.new)
    q = x %*% v.new; l.new = t.loss(q, v.new); diff = l.old - l.new

    # Display progress if verbose
    if (verbose & (i %% verbose == 0))
      progress.str(pline(i, l.new, l.old, diff))

    # Break if improvement smaller than tol, that is, sufficient convergence
    if (diff / l.old < loss.tol) break
  }

  # Ensure information of last iteration is displayed
  if (verbose & (i %% verbose)) progress.str(pline(i, l.new, l.old, diff))

  # Return mlfit.bin.fit object
  res = list('coefficients'=v.new, 'loss'=l.new, 'lambda'=lambda,
    'hinge.loss'=loss, 'huber.k'=huber.k)
  class(res) = 'mlkit.bin.fit'
  return(res)
}
