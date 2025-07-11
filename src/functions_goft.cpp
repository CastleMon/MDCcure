#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix weight_cpp(NumericVector xx, double h) {
  int n = xx.size();
  NumericMatrix w(n, n);
  double z;

  for (int i = 0; i < n; ++i) {
    double sumw = 0.0;
    for (int j = 0; j < n; ++j) {
      z = (xx[i] - xx[j]) / h;
      if (std::abs(z) > 1.0) {
        w(i, j) = 0.0;
      } else {
        w(i, j) = 0.75 * (1.0 - z * z);
      }
      sumw += w(i, j);
    }
    for (int j = 0; j < n; ++j) {
      w(i, j) /= sumw;
    }
  }
  return w;
}

// [[Rcpp::export]]
NumericVector beran_cpp(NumericMatrix data, NumericMatrix w) {
  int n = data.nrow();

  NumericMatrix k(n, n);
  for (int i = 0; i < n; i++) {
    k(i, 0) = 0.0;
    for (int j = 1; j < n; j++) {
      k(i, j) = w(i, j - 1);
    }
  }

  NumericMatrix cumsumk(n, n);
  for (int i = 0; i < n; i++) {
    double running_sum = 0.0;
    for (int j = 0; j < n; j++) {
      running_sum += k(i, j);
      cumsumk(i, j) = running_sum;
    }
  }

  NumericMatrix cumw(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cumw(i, j) = 1.0 - cumsumk(i, j);
      if (cumw(i, j) == 0.0) {
        w(i, j) = 0.0;
        cumw(i, j) = 1.0;
      }
    }
  }

  NumericMatrix delta(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      delta(i, j) = data(j, 1);
    }
  }

  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double prod = 1.0;
    for (int j = 0; j < n; j++) {
      double q = 1.0 - w(i, j) / cumw(i, j);
      if (delta(i, j) == 0.0) q = 1.0;
      prod *= q;
    }
    out[i] = prod;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix beranT_cpp(NumericMatrix data, NumericMatrix w) {
  int n = data.nrow();
  NumericMatrix k(n, n);
  NumericMatrix cumsumk(n, n);
  NumericVector cumw(n);
  NumericMatrix w_mod = clone(w);
  NumericMatrix q(n, n);

  for (int i = 0; i < n; ++i) {
    k(i, 0) = 0.0;
    for (int j = 1; j < n; ++j) {
      k(i, j) = w(i, j - 1);
    }
  }

  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int j = 0; j < n; ++j) {
      s += k(i, j);
      cumsumk(i, j) = s;
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cumw[j] = 1.0 - cumsumk(i, j);
    }

    for (int j = 0; j < n; ++j) {
      if (cumw[j] == 0.0) {
        cumw[j] = 1.0;
        w_mod(i, j) = 0.0;
      }
    }

    for (int j = 0; j < n; ++j) {
      q(i, j) = 1.0 - w_mod(i, j) / cumw[j];
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (data(j, 1) == 0.0) {
        q(i, j) = 1.0;
      }
    }
  }

  NumericMatrix out(n, n);
  for (int i = 0; i < n; ++i) {
    double prod = 1.0;
    for (int j = 0; j < n; ++j) {
      prod *= q(i, j);
      out(i, j) = prod;
    }
  }

  return out;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix beranC_cpp(NumericMatrix data, NumericMatrix w) {
  int n = data.nrow();

  NumericMatrix k(n, n);
  for (int i = 0; i < n; ++i) {
    k(i, 0) = 0.0;
    for (int j = 1; j < n; ++j) {
      k(i, j) = w(i, j - 1);
    }
  }

  NumericMatrix cumsumk(n, n);
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int j = 0; j < n; ++j) {
      s += k(i, j);
      cumsumk(i, j) = s;
    }
  }

  NumericMatrix cumw(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cumw(i, j) = 1.0 - cumsumk(i, j);
    }
  }

  NumericMatrix w_mod = clone(w);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (cumw(i, j) < 1e-10) {
        cumw(i, j) = 1.0;
        w_mod(i, j) = 0.0;
      }
    }
  }

  NumericMatrix q(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      q(i, j) = 1.0 - w_mod(i, j) / cumw(i, j);
    }
  }

  NumericMatrix delta(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      delta(i, j) = data(j, 1);
    }
  }

  NumericMatrix q_mod = clone(q);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (delta(i, j) == 1.0) {
        q_mod(i, j) = 1.0;
      }
    }
  }

  NumericMatrix out(n, n);
  for (int i = 0; i < n; ++i) {
    double prod = 1.0;
    for (int j = 0; j < n; ++j) {
      prod *= q_mod(i, j);
      out(i, j) = prod;
    }
  }

  NumericMatrix out_t(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      out_t(i, j) = out(j, i);
    }
  }

  NumericMatrix result(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result(i, j) = 1.0 - out_t(i, j);
    }
  }

  for (int j = 0; j < n; ++j) {
    result(n - 1, j) = 1.0;
  }

  return result;
}
