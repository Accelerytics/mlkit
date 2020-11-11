################################################################################
# Helper functions for computing summary statistics for linear regression
# models.
#
# Inputs:
#   beta:         Vector of parameter estimates.
#   x:            Table containing numerical explanatory variables.
#   y:            Column containing a numerical dependent variable.
#
# Output:
#   Summary statistic for given linear regression model.
rss = function(beta, x, y) {sum((y - x %*% beta) ^ 2)}
r2 = function(beta, x, y) 1 - rss(beta, x, y) / sum((y - mean(y)) ^ 2)
adj.r2 = function(beta, x, y) {
  N = nrow(x); return(1 - (1 - r2(beta, x, y)) * (N - 1) / (N - sum(beta != 0)))
}

################################################################################
# Helper functions for initializing machine learning models.
#
# Inputs:
#   beta.init:    Initial beta parameter.
#   x:            Table containing numerical explanatory variables.
#   y:            Column containing a numerical dependent variable.
#   intercept:    Indicator for whether to include an intercept or not.
#   standardize:  Indicator for whether to standardize the input. Ignored when
#                 intercept is TRUE.
#
# Output:
#   Dependent or explanatory variables formatted to be used in a regression
#   model.
create.x = function(x, intercept, standardize) {
  x = data.matrix(x)
  if (intercept) return(cbind(1, x)) else if (standardize) x = scale(x)
  return(x)
}
create.y = function(y, intercept, standardize) {
  y = data.matrix(y)
  if (!intercept & standardize) y = scale(y)
  return(y)
}
initialize.beta = function(beta.init, x) {
  if(is.null(beta.init)) return(2 * stats::runif(ncol(x)) - 1)
  return(beta.init)
}

################################################################################
# Helper functions for computing the loss of a linear regression model with an
# elastic net penalty on the parameter estimates.
#
# Inputs:
#   x:            Table containing numerical explanatory variables.
#   y:            Column containing a numerical dependent variable.
#   intercept:    Indicator for whether to include an intercept or not.
#   lambda.l1:    Penalty term for the L1-norm term.
#   lambda.l2:    Penalty term for the L2-norm term.
#
# Output:
#   Loss of the given linear model with an elastic net penalty.
elastic.net.loss = function(beta, x, y, intercept, lambda.l1, lambda.l2)
  rss(beta, x, y) / (2 * nrow(x)) +
  lambda.l1 * sum(abs(beta[(1 + intercept):ncol(x)])) +
  lambda.l2 * sum(beta[(1 + intercept):ncol(x)] ^ 2) / 2

################################################################################
# Helper function for descaling a linear regression estimator estimated on
# scaled data.
#
# Inputs:
#   beta:         Vector of parameter estimates.
#   x.mean:       Vector of column means of the explanatory variables.
#   x.sd:         Vector of column standard deviations of the explanatory
#                 variables.
#   y.mean:       Mean of dependent variable.
#   y.sd:         Standard deviation of dependent variable.
#
# Output:
#   Descaled linear regression estimator.
descale.beta = function(beta, x.mean, x.sd, y.mean, y.sd) {
  beta = y.sd * beta / x.sd
  beta = c(y.mean - sum(x.mean * beta), beta)
  return(beta)
}

################################################################################
# Helper functions for displaying progress in custom implementation of machine
# learning models.
progress.str = function(line.list) {
  for (l in line.list) cat(progress.line(l)); cat('\n\n')
}

progress.line = function(line) paste0(
  format(paste0(line[1], ':'), width=25), format(line[2], width=25,
  justify='right', nsmall=ifelse(is.integer(line[2]), 0, 10)),'\n')
