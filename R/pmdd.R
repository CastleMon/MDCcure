#' @title Partial Martingale Difference Divergence (pMDD)
#'
#' @description
#' \code{pmdd} measures conditional mean dependence of \code{Y} given \code{X}, adjusting for the dependence on \code{Z}.
#'
#' @param X A vector or matrix where rows represent samples and columns represent variables.
#' @param Y A vector or matrix where rows represent samples and columns represent variables.
#' @param Z A vector or matrix where rows represent samples and columns represent variables.
#'
#' @return Returns the squared partial martingale difference divergence of \code{Y} given \code{X}, adjusting for the dependence on \code{Z}.
#'
#' @references
#' Park, T., Shao, X., and Yao, S. (2015). Partial martingale difference correlation. \emph{Electronic Journal of Statistics}, \bold{9}(1), 1492-1517. \doi{10.1214/15-EJS1047}.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 50
#' x <- matrix(rnorm(n * 5), nrow = n)  # explanatory variables
#' y <- matrix(rnorm(n), nrow = n)      # response variable
#' z <- matrix(rnorm(n * 2), nrow = n)  # conditioning variables
#'
#' # Compute partial MDD
#' pmdd(x, y, z)
#'
#' @export
pmdd <- function(X, Y, Z) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- as.matrix(Z)
  pmdd_value <- pmdd_cpp(X, Y, Z)
  return(pmdd_value)
}
