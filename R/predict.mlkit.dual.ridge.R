#' Model Predictions for Dual Ridge Regression Model
#'
#' @description Custom implementation for prediction using model fits using the
#' \code{\link{dual.ridge.lm}} method.
#'
#' @method predict mlkit.dual.ridge.fit
#'
#' @param object \code{mlkit.dual.ridge.fit} object generated by a call to the
#' \code{\link{dual.ridge.lm}} method that is used for prediction.
#' @param newdata optional matrix of explanatory variables to use in prediction.
#' Default is \code{NULL} in which case the in-sample predictions are returned.
#' @param ... additional arguments affecting the predictions produced.
#'
#' @return Atomic vector containing predictions based on the given model and
#' explanatory variables.
#'
#' @export
#'
predict.mlkit.dual.ridge.fit = function(object, newdata=NULL, ...) {

  # Generate transformed explanatory variables matrix k
  k = create.k(object$x, object$kernel, object$const, object$degree,
    object$scale, length(object$yhat), y=newdata)

  # Generate predictions and return them in an atomic vector
  yhat = object$coefficients[1] + k %*% MASS::ginv(object$ker.mat) %*%
    object$coefficients[-1]
  return(yhat)
}
