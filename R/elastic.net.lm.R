elastic.net.lm = function(x, y, lambda, alpha, intercept=F, standardize=T,
  beta.init=NULL, beta.tol=0, loss.tol=1e-6, eps=1e-6, seed=NULL, verbose=0) {
  # Implementation of the MM algorithm solver for a linear regression model
  # an elastic net penalty term.
  #
  # Inputs:
  #   x:            Table containing numerical explanatory variables.
  #   y:            Column containing a numerical dependent variable.
  #   lambda:       Penalty scaling constant.
  #   alpha:        Scalar of penalty for L1-norm of beta. Note that the scalar
  #                 assigned to the L2-norm is equal to (1 - alpha) / 2.
  #   intercept:    Indicator for whether or not to add an intercept. Default
  #                 is FALSE. If TRUE, standardize is ignored.
  #   standardize:  Indicator for whether or not to scale data. If intercept is
  #                 TRUE, this argument is ignored. Default is TRUE.
  #   beta.init:    Optional initial value of betas. If NULL, a random beta is
  #                 sampled from the continuous uniform distribution on the
  #                 interval [0, 1]. Default is NULL.
  #   beta.tol:     Rounding tolerance for beta. Default is 0.
  #   loss.tol:     Tolerated loss, default is 1e-6.
  #   eps:          Correcting value in computation of D matrix, default is
  #                 1e-6.
  #   seed:         Optional seed, used for generating the random initial beta.
  #                 Ignored when an initial beta is provided. Default is NULL.
  #   verbose:      Integer indicating the step-size of printing progress
  #                 updates, default is 0, that is, no progress updates.
  #
  # Output:
  #   Dataframe containing the results of the linear regression model with
  #   elastic net penalty term solved using the MM algorithm.

  # Estimate model with Ridge regression if possible
  if (alpha == 0) return(ridge.lm(x, y, lambda, intercept, standardize,
    beta.tol, verbose))

  # Define descaling function
  descale = function(beta) descale.beta(beta, x, y)

  # Add intercept or standarize data if necessary
  x = create.x(x, intercept, standardize)
  y = create.y(y, intercept, standardize)

  # Define constants
  N = nrow(x); P = ncol(x); if (!is.null(seed)) set.seed(seed)
  lambda.l1 = lambda * alpha; lambda.l2 = lambda * (1 - alpha)
  lambda.l2.I = diag(rep(N * lambda.l2, P)); if (intercept) lambda.l2.I[1,1] = 0
  Xt.X = crossprod(x); Xt.y = crossprod(x, y)

  # Define helper functions for computing specific expressions and loss
  lambda.l1.D = function(beta) {
    beta[abs(beta) < eps] = eps
    lambda.l1.D = diag(N * lambda.l1 / abs(beta))
    if (intercept) lambda.l1.D[1,1] = 0
    return(lambda.l1.D)
  }
  loss = function(beta)
    elastic.net.loss(beta, x, y, intercept, lambda.l1, lambda.l2)
  pline = function(i, o, n, d)
    list(c('Iteration', i), c('Loss.old', o), c('Loss.new', n), c('Delta', d))

  # Choose some inital beta_0 and compute initial loss
  b.new = initialize.beta(beta.init, x); l.new = loss(b.new)

  # Update iteration and replace old parameters by previous until convergence
  i = 0L; while (TRUE) { i = i + 1L; l.old = l.new; b.old = b.new

    # Update parameters, loss and delta
    b.new = solve(Xt.X + lambda.l1.D(b.old) + lambda.l2.I, Xt.y)
    l.new = loss(b.new); diff = l.old - l.new

    # Display progress if verbose
    if (verbose & (i %% verbose == 0)) progress(pline(i, l.new, l.old, diff))

    # Break if improvement smaller than tol, that is, sufficient convergence
    if (diff / l.old < loss.tol) break
  }

  # Ensure information of last iteration is displayed
  if (verbose & (i %% verbose)) progress(pline(i, l.new, l.old, diff))

  # Force elements smaller than beta.tol to zero
  b.new[abs(b.new) < beta.tol] = 0

  # Descale beta if necessary
  if (intercept | !standardize) beta = c(0, b.new) else beta = descale(b.new)

  return(list('a0'=beta[1], 'beta'=beta[-1], 'alpha'=alpha, 'lambda'=lambda,
    'loss'=loss(b.new), 'R^2'=r2(b.new, x, y),
    'adjusted R^2'=adj.r2(b.new, x, y)))
}
