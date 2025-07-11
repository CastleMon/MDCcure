// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>
#include <numeric>
#include <thread>
#include <RcppParallel.h>

using namespace Rcpp;
using std::vector;
using std::string;

vector<vector<double>> numericMatrixToVector(const NumericMatrix& mat) {
  int n = mat.nrow();
  int p = mat.ncol();
  vector<vector<double>> vec(n, vector<double>(p));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < p; ++j)
      vec[i][j] = mat(i, j);
  return vec;
}

vector<double> numericVectorToVector(const NumericVector& vec) {
  int n = vec.size();
  vector<double> out(n);
  for (int i = 0; i < n; ++i)
    out[i] = vec[i];
  return out;
}

vector<vector<double>> dist_matrix_cpp(const vector<vector<double>>& X) {
  int n = X.size();
  int p = X[0].size();
  vector<vector<double>> D(n, vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < p; ++k) {
        double diff = X[i][k] - X[j][k];
        sum += diff * diff;
      }
      double d = std::sqrt(sum);
      D[i][j] = d;
      D[j][i] = d;
    }
  }
  return D;
}

vector<vector<double>> dist_vector_cpp(const vector<double>& Y) {
  int n = Y.size();
  vector<vector<double>> D(n, vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
      double d = std::abs(Y[i] - Y[j]);
      D[i][j] = d;
      D[j][i] = d;
    }
  }
  return D;
}

vector<vector<double>> u_center_cpp(const vector<vector<double>>& A) {
  int n = A.size();
  vector<double> R(n, 0.0), C(n, 0.0);
  double S = 0.0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R[i] += A[i][j];
      C[j] += A[i][j];
      S += A[i][j];
    }
  }

  double scalar = S / ((n - 1.0) * (n - 2.0));
  vector<vector<double>> UA(n, vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      UA[i][j] = A[i][j] - R[i] / (n - 2.0) - C[j] / (n - 2.0) + scalar;
    }
    UA[i][i] = 0.0;
  }
  return UA;
}

vector<vector<double>> d_center_cpp(const vector<vector<double>>& A) {
  int n = A.size();
  vector<double> R(n, 0.0), C(n, 0.0);
  double S = 0.0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R[i] += A[i][j];
      C[j] += A[i][j];
      S += A[i][j];
    }
  }

  double scalar = S / (n * n);
  vector<vector<double>> DA(n, vector<double>(n, 0.0));

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      DA[i][j] = A[i][j] - R[i] / n - C[j] / n + scalar;

  return DA;
}

double inner_u_cpp(const vector<vector<double>>& A, const vector<vector<double>>& B) {
  int n = A.size();
  double total = 0.0;

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      total += A[i][j] * B[i][j];

  return total / (n * (n - 3.0));
}

double inner_d_cpp(const vector<vector<double>>& A, const vector<vector<double>>& B) {
  int n = A.size();
  double total = 0.0;

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      total += A[i][j] * B[i][j];

  return total / (n * n);
}

double mdc_cpp_vec(const vector<vector<double>>& X, const vector<double>& Y, string center = "U") {
  int n = Y.size();

  vector<vector<double>> DX = dist_matrix_cpp(X);
  vector<vector<double>> DY = dist_vector_cpp(Y);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      DY[i][j] = 0.5 * DY[i][j] * DY[i][j];

  if (center == "U") {
    vector<vector<double>> AX = u_center_cpp(DX);
    vector<vector<double>> AY = u_center_cpp(DY);
    double num = inner_u_cpp(AX, AY);
    double den = std::sqrt(inner_u_cpp(AX, AX) * inner_u_cpp(AY, AY));
    return num / den;
  } else if (center == "D") {
    vector<vector<double>> AX = d_center_cpp(DX);
    vector<vector<double>> AY = d_center_cpp(DY);
    double num = inner_d_cpp(AX, AY);
    double den = std::sqrt(inner_d_cpp(AX, AX) * inner_d_cpp(AY, AY));
    return num / den;
  } else {
    throw std::invalid_argument("center must be 'U' or 'D'");
  }
}

struct PermutationWorkerCppVec : public RcppParallel::Worker {
  const vector<vector<double>>& X;
  const vector<double>& Y;
  const string center;
  const double observed;
  const int n;
  RcppParallel::RVector<int> results;

  std::mt19937 gen;

  PermutationWorkerCppVec(const vector<vector<double>>& X_,
                          const vector<double>& Y_,
                          const string& center_,
                          double observed_,
                          int n_,
                          Rcpp::IntegerVector& results_)
    : X(X_), Y(Y_), center(center_), observed(observed_), n(n_), results(results_) {
    std::random_device rd;
    gen.seed(rd());
  }

  void operator()(std::size_t begin, std::size_t end) {
    std::mt19937 local_gen(gen());
    vector<double> permuted_Y = Y;

    for (std::size_t i = begin; i < end; i++) {
      std::shuffle(permuted_Y.begin(), permuted_Y.end(), local_gen);
      double stat = mdc_cpp_vec(X, permuted_Y, center);
      results[i] = (stat >= observed) ? 1 : 0;
    }
  }
};

// [[Rcpp::export]]
List permutation_test_cpp_parallel(NumericMatrix X, NumericVector Y,
                                       int n_permutations,
                                       std::string center,
                                       bool parallel,
                                       int n_threads) {
  int n = Y.size();

  vector<vector<double>> X_vec = numericMatrixToVector(X);
  vector<double> Y_vec = numericVectorToVector(Y);

  double observed = mdc_cpp_vec(X_vec, Y_vec, center);

  IntegerVector results(n_permutations);

  if (parallel) {
    if (n_threads <= 0) {
      n_threads = std::max(1u, std::thread::hardware_concurrency() - 1);
    }
    PermutationWorkerCppVec worker(X_vec, Y_vec, center, observed, n_permutations, results);
    RcppParallel::parallelFor(0, n_permutations, worker, n_threads);
  } else {
    std::mt19937 gen(std::random_device{}());
    vector<double> permuted_Y = Y_vec;

    for (int i = 0; i < n_permutations; ++i) {
      std::shuffle(permuted_Y.begin(), permuted_Y.end(), gen);
      double stat = mdc_cpp_vec(X_vec, permuted_Y, center);
      results[i] = (stat >= observed) ? 1 : 0;
    }
  }

  double pvalue = (std::accumulate(results.begin(), results.end(), 0.0) + 1.0) / (n_permutations + 1.0);

  return List::create(Named("statistic") = observed,
                      Named("p.value") = pvalue);
}
