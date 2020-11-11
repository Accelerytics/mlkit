#' Grid search K-fold cross-validation
#'
#'@description
#' Implementation of the grid search approach using K-fold cross-validation for
#' hyperparameter tuning of a given \code{estimator}.
#'
#' @param formula an object of class \code{\link[stats]{formula}}: a symbolic
#' description of the model to be fitted following the standard of
#' \code{\link[stats]{lm}}.
#' @param data an optional data frame, list or environment (or object coercible
#' by \code{\link[base]{as.data.frame}} to a data frame) containing the
#' variables in the model. If not found in data, the variables are taken from
#' environment (\code{formula}), typically the environment from which this
#' function is called.
#' @param estimator estimator function that has arguments \code{formula} and
#' \code{data} and returns a list containing parameters estimates in the
#' component \code{coefficients}.
#' @param params.list list or vector containing hyperparameters and their
#' respective values to consider
#' @param n.folds optional number of folds (K) in cross-validation. Default is
#' 5.
#' @param ind.metric metric function taking in a numerical vector of
#' coefficients, a numerical matrix of explanatory variables and a numerical
#' vector of dependent variables that returns a performance metric for the
#' individual folds.
#' @param comb.metric optional function used to combine individual fold
#' metrics. Default is \code{\link[base]{mean}}.
#' @param fold.id optional vector containing numerical fold identifiers for
#' each row in the data. If \code{NULL}, \code{n.folds} is used and random fold
#' identifiers are constructed divided the observations equally over K folds.
#' If provided, \code{n.folds} is ignored. Default is \code{NULL}.
#' @param force optional boolean indicating whether or not to allow for errors
#' due to singularity when applying the estimator. If \code{TRUE}, the
#' individual metric of folds with hyperparameter combinations for which the
#' estimator is not able to estimate coefficients due to singularity, are set to
#' \code{Inf} and hence the errors are ignored. Default is \code{FALSE}.
#' @param verbose optional boolean indicating whether to show a progress bar.
#' Default is \code{FALSE}.
#' @param plot optional boolean indicating whether to generate heatmaps of
#' performance on the grid and displaying the estimated coefficients. Default
#' is \code{FALSE}.
#' @param heat.scale optional vector containing the hyperparameters and their
#' scales for the axes in the heatmaps. Default is \code{NULL}.
#' @param coef.lims optimal limits of the coefficients plot. Default is
#' \code{NULL}.
#' @param seed optimal seed to specify. Default is \code{NULL}.
#' @param ... additional arguments to be passed to the \code{estimator}
#' function.
#'
#' @return \code{grid.search.cross.validation} returns an object of
#' \code{\link{class}} \code{"gscv"}. An object of class \code{"gscv"} is a
#' list containing at least the following components:
#' \item{coefficients}{a named vector of optimal coefficients.}
#' \item{metric}{metric of the optimal hyperparameters.}
#' \item{params}{a named vector of optimal hyperparameters.}
#' @export
#'
grid.search.cross.validation = function(formula, data, estimator, params.list,
  n.folds=5, ind.metric, comb.metric=mean, fold.id=NULL, force=F, verbose=F,
  plot=F, heat.scale=NULL, coef.lims=NULL, seed=NULL, ...) {

  # Define constants
  y = data.matrix(data[, all.vars(formula)[1]]); set.seed(seed)
  x = stats::model.matrix(formula, data); N = nrow(x); opt.metric = Inf

  # Specify fold ids if not given and initialize metrics and ids vector
  if(is.null(fold.id)) fold.id = ((1:N) %% n.folds + 1)[sample(N, N)]
  else n.folds =length(unique(fold.id)); metrics = rep(NULL, n.folds);
  test.ids = matrix(nrow=n.folds, ncol=N)
  for (fold in 1:n.folds) test.ids[fold, ] = (fold.id == fold)

  # Create grid for cross validation search
  grid = expand.grid(params.list)
  n.combs = nrow(grid); metric = rep(NULL, n.combs)

  # If verbose, initialize progress bar
  if (verbose) pb = dplyr::progress_estimated(n.combs * n.folds)

  # Apply grid search
  for (i in 1:n.combs) {

    # Apply cross validation and if verbose, update progress bar
    for (fold in 1:n.folds) { if (verbose) pb$tick()$print()

      # Compute individual metric on fold
      metrics[fold] = tryCatch(
        ind.metric(do.call(estimator, c(list(formula=formula,
          data=data[!test.ids[fold, ], ]), as.list(grid[i, ]),
          list(...)))$coefficients, x[test.ids[fold, ], ], y[test.ids[fold, ]]),
        error = function(e) {
          warning(paste('Failed for', paste(names(params.list), '=', grid[i, ],
            collapse=', ')))
          if (force & grepl('.*singular.*', e$message)) return(Inf)
          stop(e)
        }
      )
    }

    # Combine performances on folds to overall performance
    metric[i] = comb.metric(metrics)
  }

  # Extract optimal hyperparameters
  opt.id = which.min(metric)
  opt.metric = metric[opt.id]; opt.params = grid[opt.id, ]

  # Estimate best beta
  opt.beta = do.call(estimator, c(list(x=x, y=y), as.list(opt.params),
    list(...)))$beta

  # Plot heatmaps and coefficients if required
  if (plot) {
    # Plot heatmaps
    combs = utils::combn(names(params.list), 2); grid$metric = metric
    for (i in 1:ncol(combs)) {
      col.x = combs[1, i]; col.y = combs[2, i]
      p = ggplot2::ggplot(data = grid, ggplot2::aes_string(x=col.x, y=col.y)) +
        ggplot2::geom_tile(ggplot2::aes(color=metric, fill=metric)) +
        ggplot2::ylab(latex2exp::TeX(paste0('$\\', col.y, '$'))) +
        ggplot2::xlab(latex2exp::TeX(paste0('$\\', col.x, '$')))
      if (!is.null(heat.scale))
        p = p + ggplot2::scale_y_continuous(trans=heat.scale[col.y]) +
        ggplot2::scale_x_continuous(trans=heat.scale[col.x])
      print(p)
    }

    # Plots beta estimates
    p = ggplot2::ggplot(data.frame(y=colnames(x), beta=as.vector(opt.beta)),
      ggplot2::aes(beta, y)) + ggplot2::geom_col() +
      ggplot2::ylab('Expl. variable') +
      ggplot2::xlab(latex2exp::TeX('$\\beta$'))
    if (!is.null(coef.lims)) p = p + ggplot2::xlim(coef.lims)
    print(p)
  }

  # Reformat optimal hyperparameters
  if (length(params.list) > 1) opt.params = c(opt.params)

  # Return gscv object
  res = list('coefficients'=opt.beta, 'metric'=opt.metric, 'params'=opt.params)
  class(res) = 'gscv'
  return(res)
}
