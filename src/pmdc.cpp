#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
double pmdc_cpp(NumericMatrix X, NumericMatrix Y, NumericMatrix Z) {
  int n = X.nrow();
  if (n != Y.nrow() || n != Z.nrow()) {
    stop("Dimensions of X, Y, and Z do not agree.");
  }

  NumericMatrix W(n, X.ncol() + Z.ncol());
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < X.ncol(); ++j) {
      W(i, j) = X(i, j);
    }
    for (int j = 0; j < Z.ncol(); ++j) {
      W(i, X.ncol() + j) = Z(i, j);
    }
  }

  NumericMatrix DW = dist_matrix(W);
  NumericMatrix D = u_center(DW);
  NumericMatrix DX = dist_matrix(W);
  NumericMatrix A = u_center(DX);

  NumericMatrix DY = dist_matrix(Y);

  NumericMatrix dist_Y = dist_matrix(Y);

  int n_dist = dist_Y.nrow();
  NumericMatrix DY2(n_dist, n_dist);
  for (int i = 0; i < n_dist; ++i) {
    for (int j = 0; j < n_dist; ++j) {
      DY2(i, j) = 0.5 * dist_Y(i, j) * dist_Y(i, j);
    }
  }

  NumericMatrix B = u_center(DY2);
  NumericMatrix DZ = dist_matrix(Z);
  NumericMatrix C = u_center(DZ);

  double beta = inner_u(B, C) / inner_u(C, C);
  NumericMatrix proj(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      proj(i, j) = B(i, j) - beta * C(i, j);
    }
  }

  double numerator = inner_u(proj, D);
  double denominator = sqrt(inner_u(proj, proj) * inner_u(D, D));

  double result = numerator / denominator;
  return result;
}
