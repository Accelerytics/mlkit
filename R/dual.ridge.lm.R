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
#' @return \code{dual.ridge.lm} returns an object of \code{\link{class}}
#' \code{mlkit.dual.ridge.fit}. An object of class \code{mlkit.dual.ridge.fit}
#' is a list containing at least the following components:
#' \item{coefficients}{a named vector of optimal coefficients.}
#' \item{alpha}{L1-weight hyperparameter in elastic net penalty term.}
#' \item{lambda}{penalty term scaling hyperparameter.}
#' \item{r2}{coefficient of determination for optimal coefficients.}
#' \item{ker.mat}{kernel matrix used in estimation.}
#' \item{kernel}{kernel transformation.}
#' \item{const}{constant parameter in the kernel transformation.}
#' \item{degree}{degree parameter in the kernel transformation.}
#' \item{scale}{scale parameter in the kernel transformation.}
#'
#' @export
#'
dual.ridge.lm = function(formula, data, lambda, intercept=F, standardize=F,
  kernel='lin', const=0, degree=NULL, scale=NULL) {

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

  # Generate transformed explanatory variables matrix k
  k = create.k(x, kernel, const, degree, scale)

  # Estimate model
  n = nrow(x)
  j = diag(n) - 1 / n
  k.til = j %*% k %*% j
  eig = eigen(k.til, symmetric = TRUE); quad.eig.vals = eig$value ^ 2
  inv.mat = diag(quad.eig.vals / (quad.eig.vals + lambda))
  w0 = mean(y)
  q.til = as.vector(eig$vectors %*% inv.mat %*% t(eig$vectors) %*% j %*% y)

  # Construct output
  coefficients = c('w0'=w0, q.til)
  y.hat = coefficients[1] + coefficients[-1]
  rss = sum((y - y.hat) ^ 2)
  r2 = 1 - rss / sum((y - mean(y)) ^ 2)

  # Return mlfit object
  res = list('coefficients'=coefficients, 'alpha'=0, 'lambda'=lambda, 'r2'=r2,
    'ker.mat'=k, 'kernel'=kernel, 'const'=const, 'degree'=degree, 'scale'=scale)
  class(res) = 'mlkit.dual.ridge.fit'
  return(res)
}
