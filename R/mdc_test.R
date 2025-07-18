#' @title MDC-Based Dependence Tests Between Multivariate Data and a Covariate
#'
#' @description
#' Computes dependence between a multivariate dataset \code{x} and a univariate covariate \code{y}
#' using different variants of the MDC (martingale difference correlation) test.
#'
#' @param x Vector or matrix  where rows represent samples, and columns represent variables.
#' @param y Covariate vector.
#' @param method Character string indicating the test to perform. One of:
#' \itemize{
#'   \item \code{"MDCU"}: U-centering permutation test.
#'   \item \code{"MDCV"}: Double-centering permutation test.
#'   \item \code{"FMDCU"}: Fast asymptotic test with U-centering.
#'   \item \code{"All"}: All of the above.
#' }
#' @param permutations Number of permutations. Defaults to 999.
#' @param parallel Logical. Whether to use parallel computing. Defaults to \code{TRUE}.
#' @param ncores Number of threads for parallel computing (used only if \code{parallel = TRUE}).
#'
#' @return A list containing the test results and p-values.
#'
#' @references
#' Shao, X., and Zhang, J. (2014). Martingale difference correlation...
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(50 * 5), nrow = 50)
#' y <- rbinom(50, 1, 0.5)
#' mdc_test(x, y, method = "FMDCU")
#'
#' @useDynLib MDCcure, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
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
