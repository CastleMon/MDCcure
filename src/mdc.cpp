#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
double mdc_cpp(SEXP X_, SEXP Y_, std::string center) {
  NumericMatrix X(X_);
  int n = X.nrow();

  NumericMatrix DX = dist_matrix(X);

  NumericMatrix DY;
  int ny;

  if (Rf_isMatrix(Y_)) {
    NumericMatrix Y(Y_);
    ny = Y.nrow();
    if (ny != n) stop("Number of rows of X and Y must match.");
    DY = dist_matrix(Y);
  } else if (Rf_isNumeric(Y_)) {
    NumericVector Y(Y_);
    ny = Y.size();
    if (ny != n) stop("Number of rows of X and Y must match.");
    DY = dist_vector(Y);
  } else {
    stop("Y must be a numeric vector or matrix.");
  }

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      DY(i, j) = 0.5 * DY(i, j) * DY(i, j);

  if (center == "U") {
    NumericMatrix AX = u_center(DX);
    NumericMatrix AY = u_center(DY);
    double num = inner_u(AX, AY);
    double den = sqrt(inner_u(AX, AX) * inner_u(AY, AY));
    return num / den;
  } else if (center == "D") {
    NumericMatrix AX = d_center(DX);
    NumericMatrix AY = d_center(DY);
    double num = inner_d(AX, AY);
    double den = sqrt(inner_d(AX, AX) * inner_d(AY, AY));
    return num / den;
  } else {
    stop("center must be 'U' or 'D'");
  }
}

