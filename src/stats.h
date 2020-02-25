#ifndef STATS_H
#define STATS_H

#include <RcppArmadillo.h>
#include <omp.h>
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <R_ext/Print.h>
#include <boost/format.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory, BH)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

NumericVector freq_hap(SEXP pBigMat, const IntegerVector index_, const int threads = 0);
NumericVector freq(SEXP pBigMat, const IntegerVector index_, const int threads = 0);
SEXP BigStat(SEXP pBigMat, const IntegerVector index_, const int threads = 0);
IntegerVector which_c(const NumericVector x, const double value, const int c = 1);
#endif
