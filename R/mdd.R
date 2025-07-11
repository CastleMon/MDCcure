#' @title Martingale Difference Divergence (MDD)
#'
#' @description
#' \code{mdd} computes the squared martingale difference divergence (MDD) between response variable(s) \code{Y}
#' and explanatory variable(s) \code{X}, measuring conditional mean dependence.
#'
#' @param X A vector or matrix where rows represent samples and columns represent variables.
#' @param Y A vector or matrix where rows represent samples and columns represent variables.
#' @param center Character string indicating the centering method to use. One of:
#' \itemize{
#'   \item \code{"U"}: U-centering, which provides an unbiased estimator.
#'   \item \code{"D"}: Double-centering, which leads to a biased estimator.
#' }
#' Default is \code{"U"}.
#'
#' @return Returns the squared Martingale Difference Divergence of \code{Y} given \code{X}.
#'
#' @references
#' Shao, X., and Zhang, J. (2014). Martingale difference correlation and its use in high-dimensional variable screening.
#'  \emph{Journal of the American Statistical Association}, \bold{109}(507), 1302-1318. \url{http://dx.doi.org/10.1080/01621459.2014.887012}.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 50
#' x <- matrix(rnorm(n * 5), nrow = n)  # multivariate explanatory variables
#' y_vec <- rbinom(n, 1, 0.5)           # univariate response
#' y_mat <- matrix(rnorm(n * 2), nrow = n)  # multivariate response
#'
#' # Compute MDD with vector Y and U-centering
#' mdd(x, y_vec, center = "U")
#'
#' # Compute MDD with matrix Y and double-centering
#' mdd(x, y_mat, center = "D")
#'
#' @export


mdd <- function(X, Y, center = "U") {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  mdd_value <- mdd_cpp(X, Y, center)
  return(mdd_value)
}
