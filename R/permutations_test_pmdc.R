permutation_test_pmdc.cpp <- function(x, y, z, n_permutations = 999) {

  set.seed(123)
  observed_statistic <- pmdc_cpp(x, y, z)

  permuted_statistics <- sapply(1:n_permutations, function(i) {
    x_perm <- x[sample(seq_len(nrow(x))), , drop = FALSE]
    pmdc_cpp(x_perm, y, z)
  })

  p_value <- (1 + sum(permuted_statistics >= observed_statistic)) / (1 + n_permutations)

  return(list(statistic = observed_statistic, p.value = p_value))
}
