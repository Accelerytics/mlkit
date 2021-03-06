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
#   kernel:       Kernel transformation function.
#   const:        Constant parameter to use in kernel transformation function.
#   degree:       Degree parameter to use in kernel transformation function.
#   scale:        Scale parameter to use in kernel transformation function.
#   n:            Number of rows parameter to use in kernel transformation
#                 function. Default is NULL in which case the number of rows of
#                 x is used.
#
# Output:
#   Dependent or explanatory variables formatted to be used in a regression
#   model.

KERNELS = c(
  'cau', # Cauchy kernel (scale)
  'cir', # Circular kernel (scale)
  'pol', # Polynomial kernel (const, scale, degree)
  'imq', # Inverse multiquadratic kernel (const)
  'lin', # Linear kernel (const)
  'log', # Logarithmic kernel (degree)
  'mqk', # Multiquadric kernel (const)
  'pow', # Power kernel (degree)
  'rbf', # Radial basis kernel (scale)
  'rqk', # Rational quadratic kernel
  'sig', # Hyperbolic tangent or sigmoid kernel (const, scale)
  'sph' # Spherical kernel (scale)
)
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
create.k = function(x, kernel, const, degree, scale, n=NULL, y=NULL) {

  # Define custom function for Euclidean distance between columns of matrices
  eucl.dist = function(x, y) outer(rowSums(y ^ 2), rep(1, nrow(x))) +
    outer(rep(1, nrow(y)), rowSums(x ^ 2)) - 2 * y %*% t(x)

  # Correct and check input
  if (is.null(y)) y = x; if (is.null(n)) n = nrow(x)
  kernel = match.arg(kernel, KERNELS, several.ok = F)

  # Return kernel transformation
  if (kernel == 'cau') {
    if (is.null(scale)) stop('Cauchy kernel requires scale argument.')
    return(1 + as.matrix(eucl.dist(x, y) ^ 2 / scale ^ 2))
  }
  if (kernel == 'cir') {
    if (is.null(scale)) stop('Circular kernel requires scale argument.')
    std.dist = as.matrix(eucl.dist(x, y) / scale)
    return((std.dist < 1) * 2 / pi * (acos(-std.dist) - std.dist * sqrt(1 -
      std.dist ^ 2)))
  }
  if (kernel == 'pol') {
    if (is.null(const)) stop('Polynomial kernel requires const argument.')
    if (is.null(degree)) stop('Polynomial kernel requires degree argument.')
    if (is.null(scale)) stop('Polynomial kernel requires scale argument.')
    return((const + scale * y %*% t(x)) ^ degree)
  }
  if (kernel == 'imq') {
    if (is.null(const))
      stop('Inverse multiquadratic kernel requires const argument.')
    return(sqrt(as.matrix(eucl.dist(x, y) ^ 2 + const ^ 2)) ^ -1)
  }
  if (kernel == 'log') {
    if (is.null(degree)) stop('Logarithmic kernel requires degree argument.')
    return(-log(as.matrix(eucl.dist(x, y) ^ degree) + 1))
  }
  if (kernel == 'mqk') {
    if (is.null(const)) stop('Multiquadratic kernel requires const argument.')
    return(sqrt(as.matrix(eucl.dist(x, y) ^ 2 + const ^ 2)))
  }
  if (kernel == 'pow') {
    if (is.null(degree)) stop('Power kernel requires degree argument.')
    return((y %*% t(x)) ^ degree)
  }
  if (kernel == 'rbf') {
    if (is.null(scale))
      stop('Radial basis function kernel requires scale argument.')
    return(exp(-scale * as.matrix(eucl.dist(x, y) ^ 2)))
  }
  if (kernel == 'rqk') {
    if (is.null(const))
      stop('Rational quadratic kernel requires const argument.')
    quad.dist = as.matrix(eucl.dist(x, y) ^ 2)
    return(1 - quad.dist / (quad.dist + const))
  }
  if (kernel == 'sig') {
    if (is.null(const)) stop('Sigmoid kernel requires const argument.')
    if (is.null(scale)) stop('Sigmoid kernel requires scale argument.')
    return(tanh(const + scale * y %*% t(x)))
  }
  if (kernel == 'sph') {
    if (is.null(scale)) stop('Spherical kernel requires scale argument.')
    std.dist = as.matrix(eucl.dist(x, y) / scale)
    return((std.dist < 1) * (1 - 3 / 2 * std.dist + std.dist ^ 3 / 2))
  }
  # Default: linear kernel
  if (is.null(const)) stop('Linear kernel requires const argument.')
  return(y %*% t(x) + const)
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
# Helper function for hinge loss function.
#
# Inputs:
#   loss:         Character hinge loss type.
#   huber.k:      Huber hinge loss hyperparameter k.
#
# Output:
#   Loss of the given linear model with an elastic net penalty.
hinge.loss = function(loss, huber.k) {
  if (loss == 'abs') return(function(x) x)
  if (loss == 'qua') return(function(x) x ^ 2)
  if (loss == 'hub') return(function(x) ifelse(abs(x) < huber.k + 1, x ^ 2 /
    (2 * (huber.k + 1)), x - (huber.k + 1) / 2))
  stop('Hinge loss type not recognized.')
}

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

################################################################################
# Helper functions for Support Vector Machines.

svm.upt.a = function(q, ind, loss, huber.k) {
  if (loss == 'abs') return(diag(as.vector(0.25 / abs(1 + (-1) ^ ind * q))))
  if (loss == 'qua') return(diag(length(ind)))
  if (loss == 'hub') return(diag(rep(0.5 / (huber.k + 1), length(ind))))
  stop('Unknown hinge loss function.')
}
svm.upt.b = function(a, q, ind, loss, huber.k) {
  if (loss == 'abs') return((-1) ^ ind * diag(a) + 0.25)
  if (loss == 'qua') return(ifelse(ind, q ^ (q > 1), ifelse(q <= -1, q, -1)))
  if (loss == 'hub') {
    if (huber.k < -1) stop("Invalid k, must be larger or equal to -1.")
    return((-1) ^ ind * diag(a) * q ^ ((q < -huber.k ^ ind) || (q > huber.k ^
      !ind)) + 0.5 * ((-1) ^ ind * q <= -huber.k))
  }
  stop('Unknown hinge loss function.')
}
