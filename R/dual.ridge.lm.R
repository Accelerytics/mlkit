#' Fitting Linear Models with Ridge Penalty using Dual Solution
#'
#' @description Implementation of the dual analytical solution for a linear
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
#' @param kernel optional kernel to use in the ridge regression model. By
#' default the linear kernel with constant zero is used, that is, no kernel
#' transformation is applied. See the details section for more details on
#' available kernel transformations.
#' @param const optional constant parameter for the kernel transformation.
#' Default is \code{0}.
#' @param degree optional degree parameter in the kernel transformation. Default
#' is \code{NULL}.
#' @param scale optional scale parameter in the kernel transformation. Default
#' is \code{NULL}.
#'
#' @return \code{ridge.lm} returns an object of \code{\link{class}}
#' \code{mlfit}. An object of class \code{mlfit} is a list containing at
#' least the following components:
#' \item{coefficients}{a named vector of optimal coefficients.}
#' \item{loss}{residual sum of squares plus ridge loss for optimal
#' coefficients.}
#' \item{r2}{coefficient of determination for optimal coefficients.}
#' @export
#'
dual.ridge.lm = function(formula, data, lambda, intercept=F, standardize=F,
  kernel='lin', const=0, degree=NULL, scale=NULL, ...) {

  # Extract dependent variable and explanatory variables
  x = stats::model.matrix(formula, data)
  y = data.matrix(data[, all.vars(formula)[1]]);

  # Store scaling parameters
  if (!intercept & !standardize) {
    x.mean = colMeans(x); x.sd = apply(x, 2, stats::sd)
    y.mean = mean(y); y.sd = stats::sd(y)
  }

  # Add intercept or standardize data if necessary
  x = create.x(x, intercept, standardize)
  y = create.y(y, intercept, standardize)

  KERNELS = c(
    'cir', # Circular kernel (scale)
    'cau', # Cauchy kernel (scale)
    'exp', # Exponential kernel (scale)
    'gau', # Gaussian kernel (scale)
    'pol', # Polynomial kernel (const, scale, degree)
    'imq', # Inverse multiquadratic kernel (const)
    'lap', # Laplacian kernel (scale)
    'lin', # Linear kernel (const)
    'log', # Logarithmic kernel (degree)
    'mqk', # Multiquadric kernel (const)
    'pow', # Power kernel (degree)
    'rbf', # Radial basis kernel (scale)
    'rqk', # Rational quadratic kernel
    'sig', # Hyperbolic tangent or sigmoid kernel (const, scale)
    'sph' # Spherical kernel (scale)
  )

  kernel = match.arg(kernel, KERNELS, several.ok = F)
  if (kernel == 'cau') {
    if (is.null(scale)) stop('Cauchy kernel requires scale argument.')
    k = 1 + as.matrix(stats::dist(x) ^ 2 / scale ^ 2)
  } else if (kernel == 'cir') {
    if (is.null(scale)) stop('Circular kernel requires scale argument.')
    std.dist = as.matrix(stats::dist(x) / scale)
    k = (std.dist < 1) * 2 / pi * (acos(-std.dist) - std.dist *
                                     sqrt(1 - std.dist ^ 2))
  } else if (kernel == 'exp') {
    if (is.null(scale)) stop('Exponential kernel requires scale argument.')
    k = exp(-as.matrix(stats::dist(x) / (nrow(x) * scale ^ 2)))
  } else if (kernel == 'gau') {
    if (is.null(scale)) stop('Gaussian kernel requires scale argument.')
    k = exp(-as.matrix(stats::dist(x) ^ 2 / (scale * nrow(x))))
  } else if (kernel == 'pol') {
    if (is.null(const)) stop('Polynomial kernel requires const argument.')
    if (is.null(degree)) stop('Polynomial kernel requires degree argument.')
    if (is.null(scale)) stop('Polynomial kernel requires scale argument.')
    k = (const + scale * x %*% t(x)) ^ degree
  } else if (kernel == 'imq') {
    if (is.null(const))
      stop('Inverse multiquadratic kernel requires const argument.')
    k = sqrt(as.matrix(stats::dist(x) ^ 2 + const ^ 2)) ^ -1
  } else if (kernel == 'lap') {
    if (is.null(scale)) stop('Laplacian kernel requires scale argument.')
    k = exp(as.matrix(stats::dist(x, method='manhattan') / scale))
  } else if (kernel == 'log') {
    if (is.null(degree)) stop('Logarithmic kernel requires degree argument.')
    k = -log(as.matrix(stats::dist(x) ^ degree) + 1)
  } else if (kernel == 'mqk') {
    if (is.null(const)) stop('Multiquadratic kernel requires const argument.')
    k = sqrt(as.matrix(stats::dist(x) ^ 2 + const ^ 2))
  } else if (kernel == 'pow') {
    if (is.null(degree)) stop('Power kernel requires degree argument.')
    k = (x %*% t(x)) ^ degree
  } else if (kernel == 'rbf') {
    if (is.null(scale))
      stop('Radial basis function kernel requires scale argument.')
    k = exp(-scale * as.matrix(stats::dist(x) ^ 2))
  } else if (kernel == 'rqk') {
    if (is.null(const))
      stop('Rational quadratic kernel requires const argument.')
    quad.dist = as.matrix(stats::dist(x) ^ 2)
    k = 1 - quad.dist / (quad.dist + const)
  } else if (kernel == 'sig') {
    if (is.null(const)) stop('Sigmoid kernel requires const argument.')
    if (is.null(scale)) stop('Sigmoid kernel requires scale argument.')
    k = tanh(const + scale * x %*% t(x))
  } else if (kernel == 'sph') {
    if (is.null(scale)) stop('Spherical kernel requires scale argument.')
    std.dist = as.matrix(stats::dist(x) / scale)
    k = (std.dist < 1) * (1 - 3 / 2 * std.dist + std.dist ^ 3 / 2)
  } else { # Default: linear kernel
    if (is.null(const)) stop('Linear kernel requires const argument.')
    k = x %*% t(x) + const
  }

  # Estimate model
  n = nrow(x)
  j = diag(n) - 1 / n
  k.til = j %*% k %*% j
  eig = eigen(k.til, symmetric = TRUE); quad.eig.vals = eig$value ^ 2
  inv.mat = quad.eig.vals / (quad.eig.vals + lambda)
  w0 = mean(y)
  q.til = eig$vectors %*% inv.mat %*% tcrossprod(eig$vectors, j %*% y)

  # Construct output
  coefficients = c('w0'=w0, q.til)
  y.hat = coeffients[1] + coefficients[-1]
  rss = sum((y - y.hat) ^ 2)
  r2 = 1 - rss / sum((y - mean(y)) ^ 2)

  # Return mlfit object
  res = list('coefficients'=coefficients, 'alpha'=0, 'lambda'=lambda,
    'loss'=rss + lambda * crossprod(q.til, solve(k)) %*% q.til, 'r2'=r2)
  class(res) = 'mlfit'
  return(res)
}
