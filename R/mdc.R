#' @title Martingale Difference Correlation (MDC)
#'
#' @description
#' \code{mdc} computes the squared martingale difference correlation between a response variable \code{Y}
#' and explanatory variable(s) \code{X}, measuring conditional mean dependence.
#' \code{X} can be either univariate or multivariate.
#'
#' @param X A vector or matrix where rows represent samples and columns represent variables.
#' @param Y A vector or matrix where rows represent samples and columns represent variables.
#' @param center Character string indicating the centering method to use. One of:
#' \itemize{
#'   \item \code{"U"}: U-centering, which provides an unbiased estimator.
#'   \item \code{"D"}: Double-centering, which leads to a biased estimator.
#' }
#'
#' @return Returns the squared martingale difference correlation of \code{Y} given \code{X}.
#'
#' @references
#' Shao, X., and Zhang, J. (2014). Martingale difference correlation and its use in high-dimensional variable screening.
#' \emph{Journal of the American Statistical Association}, \bold{109}(507), 1302-1318. \url{10.1080/01621459.2014.887012}.
#'
#' @seealso \code{\link{mdd}}, \code{\link{mdc_test}}
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 50
#' x <- matrix(rnorm(n * 5), nrow = n)  # multivariate data with 5 variables
#' y <- rbinom(n, 1, 0.5)               # binary covariate
#'
#' # Compute MDC with U-centering
#' mdc(x, y, center = "U")
#'
#' # Compute MDC with double-centering
#' mdc(x, y, center = "D")
#'
#' @export

mdc <- function(X, Y, center = "U"){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  mdc_value <- mdc_cpp(X, Y, center)
  return(mdc_value)
}
