#' @title MDC-Based Dependence Tests Between Multivariate Data and a Covariate
#'
#' @description
#' Computes dependence between a multivariate dataset \code{x} and a univariate covariate \code{y}
#' using different variants of the MDC (martingale difference correlation) test.
#' Available methods include permutation-based tests with U-centering (\code{"MDCU"}) or
#' double-centering (\code{"MDCV"}), a fast asymptotic version (\code{"FMDCU"}), or all of them at once (\code{"All"}).
#'
#' The permutation tests estimate p-values empirically, while the fast test computes an approximate
#' p-value based on the chi-square distribution.
#'
#' @param x Vector or matrix  where rows represent samples, and columns represent variables.
#' @param y Covariate vector.
#' @param method Character string indicating the test to perform.
#' One of:
#' \itemize{
#'   \item \code{"MDCU"}: Performs the MDC test with U-centering using permutations.
#'   \item \code{"MDCV"}: Performs the MDC test with double-centering using permutations.
#'   \item \code{"FMDCU"}: Performs the fast MDC test with U-centering.
#'   \item \code{"All"}: Performs all the above tests.
#' }
#' @param permutations Number of permutations to use for the \code{"MDCU"} and \code{"MDCV"} methods.
#' Defaults to 999.
#' @param parallel Logical. If \code{TRUE}, parallel computing is used. Defaults to \code{TRUE}.
#' @param ncores Integer. Number of threads to use for parallel computing
#' (used only if \code{parallel = TRUE}).
#' By default, it uses one less than the number of available CPU cores,
#' ensuring at least one thread is used.
#'
#' @return A list containing the results of the selected test(s), including the corresponding p-value(s).
#'
#' @references
#' Shao, X., and Zhang, J. (2014). Martingale difference correlation and its use in high-dimensional variable screening.
#' \emph{Journal of the American Statistical Association}, \bold{109}(507), 1302-1318. \doi{10.1080/01621459.2014.887012}.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 50
#' x <- matrix(rnorm(n * 5), nrow = n)  # multivariate data with 5 variables
#' y <- rbinom(n, 1, 0.5)               # binary covariate
#'
#' # Run the fast MDC test with U-centering
#' mdc_test(x, y, method = "FMDCU")
#'
#' # Run permutation-based MDC test with U-centering
#' mdc_test(x, y, method = "MDCU", permutations = 999)
#'
#' # Run all tests in parallel
#' mdc_test(x, y, method = "All", permutations = 499, parallel = TRUE)
#' @export


mdc_test <- function(x, y, method, permutations = 999, parallel = TRUE, ncores = -1){
  x <- as.matrix(x)
  result <- switch(
    method,
    "All" = {
      permMDC_U <- permutation_test_cpp_parallel(x, y, permutations, "U", parallel, ncores)

      permMDC_V <- permutation_test_cpp_parallel(x, y, permutations, "D", parallel, ncores)

      FMDCU_stat <- mdc_cpp(x, y, center = "U")
      FMDCU <- list(
        statistic = FMDCU_stat,
        p.value = 1 - pchisq(length(y) * FMDCU_stat + 1, df = 1)
      )

      list(MDCU = permMDC_U, MDCV = permMDC_V, FMDCU = FMDCU)
    },

    "MDCU" = permutation_test_cpp_parallel(x, y, permutations, "U", parallel, ncores),

    "MDCV" = permutation_test_cpp_parallel(x, y, permutations, "D", parallel, ncores),

    "FMDCU" = {
      FMDCU_stat <- mdc_cpp(x, y, center = "U")
      list(statistic = FMDCU_stat, p.value = 1 - pchisq(length(y) * FMDCU_stat + 1, df = 1))
    },
    stop("Unknown method")
  )

  # ======================= Results ==============================
  #===============================================================

  cat("============== Covariate hypothesis test for the cure rate ============\n")
  cat("Method:", method, "\n\n")

  if (method %in% c("MDCU", "MDCV")) {
    desc <- ifelse(method == "MDCU",
                   "Martingale difference correlation unbiased with permutations",
                   "Martingale difference correlation biased with permutations")
    cat(desc, "\n\n")
    cat("Number of permutations:", permutations, "\n\n")
    cat("Test statistic:", round(result$statistic, 5), "\n\n")
    cat("p-value:", signif(result$p.value, 5), "\n\n")

  } else if (method == "FMDCU") {
    cat("Fast Chi-square test based on MDC unbiased\n\n")
    cat("Test statistic:", round(result$statistic, 5), "\n\n")
    cat("p-value:", signif(result$p.value, 5), "\n\n")

  }
    else if (method == "All") {
    cat("MDCU: Martingale difference correlation unbiased with permutations\n")
    cat("Test statistic:", round(result$MDCU$statistic, 5), "\n")
    cat("p-value:", signif(result$MDCU$p.value, 5), "\n")
    cat("Number of permutations:", permutations, "\n\n")

    cat("MDCV: Martingale difference correlation biased with permutations\n")
    cat("Test statistic:", round(result$MDCV$statistic, 5), "\n")
    cat("p-value:", signif(result$MDCV$p.value, 5), "\n")
    cat("Number of permutations:", permutations, "\n\n")

    cat("FMDCU: Fast Chi-square test based on MDC unbiased\n")
    cat("Test statistic:", round(result$FMDCU$statistic, 5), "\n")
    cat("p-value:", signif(result$FMDCU$p.value, 5), "\n\n")
  }

  cat("========================================================================\n")

  invisible(list(result = result))
}
