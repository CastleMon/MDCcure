#' @export
print.testcov2 <- function(x, ...) {
  X <- x$varnames$x
  Z <- x$varnames$z

  cat("---------------------------------------------------------------------\n")
  cat("Covariate hypothesis test for the cure rate with two covariates \n")
  cat("---------------------------------------------------------------------\n")
  cat("Hypotheses:\n")
  cat(paste0("  H0: The conditional mean of the cure status given ",X, " adjusting on ", Z,
             " does not depend on ",X, ",\n"))
  cat(paste0("      i.e., E[nu|",X, ",",Z, "] = E[nu|",Z,"].\n"))
  cat(paste0("  H1: The conditional mean of the cure status depends on ", X, " adjusting on ", Z, "\n\n"))

  cat("Test Statistic:", format(x$test$statistic, digits = 4), "\n")
  cat("p-value:", format.pval(x$test$p.value, digits = 4), "\n")
  cat("Number of permutations:", x$test$B, "\n")

  invisible(x)
}
