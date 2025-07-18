// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec kernel_multivariate_cpp(const arma::rowvec& x, const arma::mat& X, const arma::mat& H) {
  int d = X.n_cols;

  double det_H = arma::det(H);
  if (det_H <= 0) {
    stop("The H matrix must be positive definite (det > 0).");
  }

  arma::mat diff = X.each_row() - x;

  arma::mat solved = diff;
  solved.each_row() /= H.diag().t();
  solved = solved.t();

  arma::vec exponent = arma::sum(diff.t() % solved, 0).t();

  double norm_const = std::pow(2 * M_PI, -0.5 * d) * std::pow(det_H, -0.5);
  arma::vec kernel_values = norm_const * arma::exp(-0.5 * exponent);

  return kernel_values;
}

arma::vec beran_estimator_multivariate_cpp(const arma::vec& t,
                                           const arma::rowvec& x,
                                           const arma::vec& T,
                                           const arma::mat& X,
                                           const arma::vec& delta,
                                           Nullable<arma::mat> H_in = R_NilValue) {
  int n = X.n_rows;
  int d = X.n_cols;

  arma::mat H;
  if (H_in.isNotNull()) {
    H = as<arma::mat>(H_in);
  } else {
    arma::vec range_d1 = {max(X.col(0)) - min(X.col(0))};
    arma::vec range_d2 = {max(X.col(1)) - min(X.col(1))};
    double h1 = 2 * range_d1(0) * std::pow(n, -1.0 / 5.0);
    double h2 = 2 * range_d2(0) * std::pow(n, -1.0 / 5.0);
    H = arma::diagmat(arma::vec({h1, h2}));
  }

  arma::vec Kh = kernel_multivariate_cpp(x, X, H);
  arma::vec Bh = Kh / arma::sum(Kh);

  // Cumulative sum in reverse
  arma::vec cumulative_Bh(n);
  double cum_sum = 0.0;
  for (int i = n - 1; i >= 0; --i) {
    cum_sum += Bh(i);
    cumulative_Bh(i) = cum_sum;
  }

  int m = t.n_elem;
  arma::vec survival(m, arma::fill::ones);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (T(j) <= t(i)) {
        double denom = cumulative_Bh(j);
        if (denom > 0) {
          survival(i) *= (1.0 - delta(j) * Bh(j) / denom);
        }
      }
    }
  }

  return survival;
}

double p_hat_cpp(const arma::rowvec& x,
                 const arma::vec& T,
                 const arma::mat& X,
                 const arma::vec& delta,
                 Nullable<arma::mat> H_in = R_NilValue) {

  arma::vec T_event = T.elem(arma::find(delta == 1));
  if (T_event.n_elem == 0) {
    stop("No events (delta == 1) in the data.");
  }
  double T1_max = T_event.max();

  arma::vec beran_value = beran_estimator_multivariate_cpp(arma::vec{T1_max}, x, T, X, delta, H_in);
  return 1.0 - beran_value(0);
}

// [[Rcpp::export]]
arma::vec latency_estimator_multivariate_cpp(const arma::vec& t,
                                             const arma::rowvec& x,
                                             const arma::vec& T,
                                             const arma::mat& X,
                                             const arma::vec& delta,
                                             Nullable<arma::mat> H_in = R_NilValue) {

  double p_hat_h = p_hat_cpp(x, T, X, delta, H_in);

  if (p_hat_h <= 0.0) {
    stop("p_hat is zero or negative, which would make the latency estimator undefined.");
  }

  arma::vec S_hat = beran_estimator_multivariate_cpp(t, x, T, X, delta, H_in);

  arma::vec S_0_b = (S_hat - (1.0 - p_hat_h)) / p_hat_h;

  return S_0_b;
}
