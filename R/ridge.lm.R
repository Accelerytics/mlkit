ridge.lm = function(x, y, lambda, intercept=F, standardize=T, descale=T,
  beta.tol=0, verbose=0) {
  # Implementation of the analytical solution for the linear regression model
  # with a Ridge penalty term.
  #
  # Inputs:
  #   x:            Table containing numerical explanatory variables.
  #   y:            Column containing a numerical dependent variable.
  #   lambda:       Penalty scaling constant.
  #   intercept:    Indicator for whether or not to add an intercept. Default
  #                 is FALSE. If TRUE, standardize is ignored.
  #   standardize:  Indicator for whether or not to scale data. Default is TRUE.
  #   beta.tol:     Rounding tolerance for beta, default is 0.
  #   verbose:      Integer indicating the step-size of printing progress
  #                 updates, default is 0, that is, no progress updates.
  #
  # Output:
  #   Dataframe containing the results of the linear regression model with
  #   Ridge penalty term solved using the analytical solution

  # Define descaling function
  descale = function(beta) descale.beta(beta, x, y)

  # Add intercept or standarize data if necessary
  x = create.x(x, intercept, standardize)
  y = create.y(y, intercept, standardize)

  # Derive beta estimate
  b.new = solve(crossprod(x) + nrow(x) * lambda / 2 * diag(ncol(x)),
    crossprod(x, y))

  # Set elements smaller than beta.tol to zero
  b.new[abs(b.new) < beta.tol] = 0

  # Descale beta if necessary
  if (intercept | !standardize) beta = c(0, b.new) else beta = descale(b.new)

  return(list('a0'=beta[1], 'beta'=beta[-1], 'alpha'=0, 'lambda'=lambda,
    'loss'=elastic.net.loss(b.new, x, y, intercept, lambda.l1=0, lambda),
    'R^2'=r2(b.new, x, y), 'adjusted R^2'=adj.r2(b.new, x, y)))
}
