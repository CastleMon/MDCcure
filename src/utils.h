#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

// Aquí defines todas las funciones auxiliares
NumericMatrix dist_matrix(NumericMatrix X);
NumericMatrix dist_vector(NumericVector x);
NumericMatrix u_center(NumericMatrix A);
NumericMatrix d_center(NumericMatrix A);
double inner_u(NumericMatrix A, NumericMatrix B);
double inner_d(NumericMatrix A, NumericMatrix B);

#endif
