#include "utils.h"

NumericMatrix dist_matrix(NumericMatrix X) {
  int n = X.nrow();
  NumericMatrix D(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < X.ncol(); ++k) {
        double diff = X(i, k) - X(j, k);
        sum += diff * diff;
      }
      double d = sqrt(sum);
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

NumericMatrix dist_vector(NumericVector x) {
  int n = x.size();
  NumericMatrix D(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double d = std::abs(x[i] - x[j]);
      D(i, j) = d;
      D(j, i) = d;
    }
  }
  return D;
}

NumericMatrix u_center(NumericMatrix A) {
  int n = A.nrow();
  NumericVector R = rowSums(A);
  NumericVector C = colSums(A);
  double S = sum(A);

  NumericMatrix UA(n, n);
  double scalar = S / ((n - 1.0) * (n - 2.0));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      UA(i, j) = A(i, j) - R[i] / (n - 2.0) - C[j] / (n - 2.0) + scalar;
    }
    UA(i, i) = 0;
  }
  return UA;
}


NumericMatrix d_center(NumericMatrix A) {
  int n = A.nrow();
  NumericVector R = rowSums(A);
  NumericVector C = colSums(A);
  double S = sum(A);

  NumericMatrix DA(n, n);
  double scalar = S / (n * n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      DA(i, j) = A(i, j) - R[i] / n - C[j] / n + scalar;
    }
  }
  return DA;
}

double inner_u(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  double total = 0.0;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      total += A(i, j) * B(i, j);

  return total / (n * (n - 3.0));
}

double inner_d(NumericMatrix A, NumericMatrix B) {
  int n = A.nrow();
  double total = 0.0;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      total += A(i, j) * B(i, j);

  return total / (n * n);
}
