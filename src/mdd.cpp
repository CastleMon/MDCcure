#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
double mdd_cpp(NumericMatrix X, NumericMatrix Y, std::string center) {
  int n = X.nrow();
  if (Y.nrow() != n) {
    stop("Dimensions of X and Y do not agree.");
  }

  NumericMatrix A = dist_matrix(X);
  NumericMatrix B = dist_matrix(Y);

  NumericMatrix B_squared(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      B_squared(i, j) = 0.5 * B(i, j) * B(i, j);
    }
  }

  if (center == "U") {
    NumericMatrix UA = u_center(A);
    NumericMatrix UB = u_center(B_squared);
    return inner_u(UA, UB);
  } else if (center == "D") {
    NumericMatrix DA = d_center(A);
    NumericMatrix DB = d_center(B_squared);
    return inner_d(DA, DB);
  } else {
    stop("Invalid center option. Use 'U' or 'D'.");
  }
}
