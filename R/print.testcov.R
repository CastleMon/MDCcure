
#' @export
print.testcov <- function(x, ...) {
  cat("============== Covariate hypothesis test for the cure rate ============\n")
  cat("Method:", x$method, "\n\n")

  result <- x$result
  method <- x$method
  P <- x$P

  if (method %in% c("MDCU", "MDCV")) {
    desc <- ifelse(method == "MDCU",
                   "Martingale difference correlation unbiased with permutations",
                   "Martingale difference correlation biased with permutations")
    cat(desc, "\n\n")
    cat("Number of permutations:", P, "\n\n")
    cat("Test statistic:", round(result$statistic, 5), "\n\n")
    cat("p-value:", signif(result$p.value, 5), "\n\n")

  } else if (method == "FMDCU") {
    cat("Fast Chi-square test based on MDC unbiased\n\n")
    cat("Test statistic:", round(result$statistic, 5), "\n\n")
    cat("p-value:", signif(result$p.value, 5), "\n\n")

  } else if (method == "GOFT") {
    cat("Number of bootstrap replicates:", P, "\n\n")
    cat("Test statistic:", round(result$statistic, 5), "\n\n")
    cat("p-value:", signif(result$p.value, 5), "\n\n")

  } else if (method == "All") {
    cat("MDCU: Martingale difference correlation unbiased with permutations\n")
    cat("Test statistic:", round(result$MDCU$statistic, 5), "\n")
    cat("p-value:", signif(result$MDCU$p.value, 5), "\n")
    cat("Number of permutations:", P, "\n\n")

    cat("MDCV: Martingale difference correlation biased with permutations\n")
    cat("Test statistic:", round(result$MDCV$statistic, 5), "\n")
    cat("p-value:", signif(result$MDCV$p.value, 5), "\n")
    cat("Number of permutations:", P, "\n\n")

    cat("FMDCU: Fast Chi-square test based on MDC unbiased\n")
    cat("Test statistic:", round(result$FMDCU$statistic, 5), "\n")
    cat("p-value:", signif(result$FMDCU$p.value, 5), "\n\n")

    cat("GOFT \n")
    cat("Test statistic:", round(result$GOFT$statistic, 5), "\n")
    cat("p-value:", signif(result$GOFT$p.value, 5), "\n")
    cat("Number of bootstrap replicates:", P, "\n\n")
  }

  cat("========================================================================\n")
  invisible(x)
}
