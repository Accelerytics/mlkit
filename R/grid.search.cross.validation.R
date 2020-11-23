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
#' predictions for the dependent variable together with a numerical vector of
#' actual outcomes of the dependent variable that returns a performance metric
#' for the individual folds.
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
#' @param contour.scale optional vector containing the hyperparameters and
#' their scales for the axes in the contour plots. Default is \code{NULL}.
#' @param coef.lims optimal limits of the coefficients plot. Default is
#' \code{NULL}.
#' @param coef.names optional names of coefficients, used in estimates barplot.
#' Default is \code{NULL} in which case the column names of the explanatory
#' variables are used.
#' @param seed optimal seed to specify. Default is \code{NULL}.
#' @param use.formula whether or not to use the formula data combination of X
#' and y as inputs to the model. Default is \code{TRUE}.
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
  plot=F, contour.scale=NULL, coef.lims=NULL, coef.names=NULL, seed=NULL,
  use.formula=T, ...) {

  # Define constants
  y = data.matrix(data[, all.vars(formula)[1]]); set.seed(seed)
  x = stats::model.matrix(formula, data); N = nrow(x); opt.metric = Inf;

  # Specify fold ids if not given and initialize metrics and ids vector
  if(is.null(fold.id)) fold.id = ((1:N) %% n.folds + 1)[sample(N, N)]
  else n.folds =length(unique(fold.id)); metrics = rep(NULL, n.folds);
  test.ids = matrix(nrow=n.folds, ncol=N)
  for (fold in 1:n.folds) test.ids[fold, ] = (fold.id == fold)

  # Create grid for cross validation search
  grid = expand.grid(params.list, stringsAsFactors=F); n.combs = nrow(grid)
  metric = rep(NULL, n.combs); n.pars = length(params.list)

  # If verbose, initialize progress bar
  if (verbose) pb = progress::progress_bar$new(
    format="[:bar] :current/:total :percent eta: :eta", total=n.combs * n.folds)

  # Apply grid search
  for (i in 1:n.combs) {

    # Specify hyperparameter combination
    params = grid[i, ]
    for (j in 1:n.pars)
      if (typeof(params[[j]]) == 'list') params[[j]] = params[[j]][[1]]

    # Apply cross validation and if verbose, update progress bar
    for (fold in 1:n.folds) { if (verbose) pb$tick()

      # Specify dependent and explanatory variables for fold
      if (use.formula)
        fold.vars = list(formula=formula, data=data[!test.ids[fold, ], ])
      else
        fold.vars = list(X=x[!test.ids[fold, ], ], y=y[!test.ids[fold, ]])

      # Compute individual metric on fold
      metrics[fold] = tryCatch(ind.metric(stats::predict(do.call(estimator,
          c(fold.vars, params, list(...))), x[test.ids[fold, ], ]),
          y[test.ids[fold, ]]),
        error = function(e) {
          warning(paste('Failed for', paste(names(params.list), '=', params,
            collapse=', ')))
          if (force & grepl('.*(singular)|(infinite).*', e$message)) return(Inf)
          stop(e)
        }
      )
    }

    # Combine performances on folds to overall performance
    metric[i] = comb.metric(metrics)
  }

  # Extract optimal hyperparameters
  opt.id = which.min(metric); opt.metric = metric[opt.id]
  opt.pars = grid[opt.id, ]
  for (j in 1:n.pars)
    if (typeof(opt.pars[[j]]) == 'list')
      opt.pars[[j]] = opt.pars[[j]][[1]]

  # Update grid to include all metrics
  grid$metric = metric

  # Estimate best beta
  if (use.formula) opt.beta = do.call(estimator, c(list(formula=formula,
    data=data), opt.pars, list(...)))$coefficients
  else opt.beta = do.call(estimator, c(list(X=x, y=y), opt.pars,
    list(...)))$coefficients

  # Plot heatmaps and coefficients if required
  if (plot) {

    # Contour plots if more than 1 hyperparameter
    if (length(params.list) > 1) {
      combs = utils::combn(names(params.list), 2);
      for (i in 1:ncol(combs)) {
        col.x = combs[1, i]; col.y = combs[2, i]
        p = ggplot2::ggplot(data = grid, ggplot2::aes_string(x=col.x, y=col.y,
          z='metric')) + ggplot2::geom_contour_filled() +
          ggplot2::ylab(latex2exp::TeX(paste0('$\\', col.y, '$'))) +
          ggplot2::xlab(latex2exp::TeX(paste0('$\\', col.x, '$')))
        if (!is.null(contour.scale))
          p = p + ggplot2::scale_y_continuous(trans=contour.scale[col.y]) +
          ggplot2::scale_x_continuous(trans=contour.scale[col.x])
      }
    } else {
      col.x = names(params.list)[1]; p = ggplot2::ggplot(data = grid,
        ggplot2::aes_string(x=col.x, y=metric)) + ggplot2::geom_line() +
        ggplot2::ylab('metric') +
        ggplot2::xlab(latex2exp::TeX(paste0('$\\', col.x, '$')))
      if (!is.null(contour.scale))
        p = p + ggplot2::scale_x_continuous(trans=contour.scale[col.x])
    }
    print(p)

    # Plots beta estimates
    if (is.null(coef.names)) coef.names = c('(Intercept)', colnames(x))
    p = ggplot2::ggplot(data.frame(y=coef.names,
      beta=as.vector(opt.beta)), ggplot2::aes(beta, y)) + ggplot2::geom_col() +
      ggplot2::ylab('Expl. variable') +
      ggplot2::xlab(latex2exp::TeX('$\\beta$'))
    if (!is.null(coef.lims)) p = p + ggplot2::xlim(coef.lims)
    print(p)
  }

  # Return gscv object
  res = list('coefficients'=opt.beta, 'metric'=opt.metric,
    'params'=c(opt.pars), 'grid'=grid)
  class(res) = 'gscv'
  return(res)
}
