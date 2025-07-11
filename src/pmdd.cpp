#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
double pmdd_cpp(NumericMatrix X, NumericMatrix Y, NumericMatrix Z) {
  int n = X.nrow();

  if (Y.nrow() != n || Z.nrow() != n)
    stop("The dimensions of X, Y, Z do not agree.");

  int p = X.ncol(), r = Z.ncol();
  NumericMatrix W(n, p + r);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      W(i, j) = X(i, j);
    }
    for (int j = 0; j < r; ++j) {
      W(i, p + j) = Z(i, j);
    }
  }

  NumericMatrix D = u_center(dist_matrix(W));
  NumericMatrix B = dist_matrix(Y);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      B(i, j) = 0.5 * B(i, j) * B(i, j);

  B = u_center(B);
  NumericMatrix C = u_center(dist_matrix(Z));

  double beta_num = inner_u(B, C);
  double beta_den = inner_u(C, C);
  double beta = beta_num / beta_den;

  // B - beta * C
  NumericMatrix proj(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      proj(i, j) = B(i, j) - beta * C(i, j);

  double result = inner_u(proj, D);
  return result;
}
