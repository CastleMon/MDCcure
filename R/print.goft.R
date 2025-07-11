#' @export
print.goft <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Goodness-of-fit test for the cure rate in a mixture cure model\n")
  cat("-------------------------------------------------------------------\n\n")
  cat("Model:", x$model, "\n\n")

  cat("Hypotheses:\n")
  cat("  H0: The model fits the data, i.e., p(x) = p_theta(x)\n")
  cat("  H1: The model does not fit the data, i.e., p(x) != p_theta(x)\n\n")

  cat("Test Statistic:", format(x$statistic, digits = 4), "\n")
  cat("p-value:", format.pval(x$p.value, digits = 4), "\n")
  cat("Number of bootstrap replications:", x$B, "\n\n")

  invisible(x)
}
