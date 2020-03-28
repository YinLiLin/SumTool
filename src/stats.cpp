// #if !defined(ARMA_64BIT_WORD)
// #define ARMA_64BIT_WORD 1
// #endif

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

/*
// [[Rcpp::export]]
arma::uvec std_setdiff(arma::uvec x, arma::uvec y){

	std::vector<int> a = arma::conv_to<std::vector<int> >::from(arma::sort(x));
	std::vector<int> b = arma::conv_to<std::vector<int> >::from(arma::sort(y));
	std::vector<int> out;

	std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(out, out.end()));

	return arma::conv_to<arma::uvec>::from(out);
}
*/


template <typename T>
bool hasNA(XPtr<BigMatrix> pMat, double NA_C, const int threads=0) {

	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    bool HasNA = false;

    #pragma omp parallel for schedule(dynamic) shared(HasNA)
    for (size_t j = 0; j < pMat->ncol(); j++) {
    	if(HasNA)	continue;
        for (size_t i = 0; i < pMat->nrow(); i++) {
            if (mat[j][i] == NA_C) {
                HasNA = true;
            }
        }
    }
    return HasNA;
}

// [[Rcpp::export]]
bool hasNA(SEXP pBigMat, const int threads=0) {
    XPtr<BigMatrix> xpMat(pBigMat);

    switch(xpMat->matrix_type()) {
    case 1:
        return hasNA<char>(xpMat, NA_CHAR, threads);
    case 2:
        return hasNA<short>(xpMat, NA_SHORT, threads);
    case 4:
        return hasNA<int>(xpMat, NA_INTEGER, threads);
    case 8:
        return hasNA<double>(xpMat, NA_REAL, threads);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

template <typename T>
NumericVector freq_hap(XPtr<BigMatrix> pMat, const IntegerVector index_, const int threads = 0){

    if (threads == 0) {
        omp_set_num_threads(omp_get_num_procs());
    }else if(threads > 0) {
        omp_set_num_threads(threads);
    }

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();
	int j, k, m = index_.size();
	double p1 = 0.0;
	NumericVector freq(m);

	#pragma omp parallel for schedule(dynamic) private(p1, k)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			if(bigm[index_[j]][k] == 12){
				p1 += 1;
			}else if(bigm[index_[j]][k] == 21){
				p1 += 1;
			}else if(bigm[index_[j]][k] == 11){
				p1 += 2;
			}
		}
		freq[j] = p1 / (ind * 2);
	}
	return freq;
}

// [[Rcpp::export]]
NumericVector freq_hap(SEXP pBigMat, const IntegerVector index_, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return freq_hap<char>(xpMat, index_, threads);
	case 2:
		return freq_hap<short>(xpMat, index_, threads);
	case 4:
		return freq_hap<int>(xpMat, index_, threads);
	case 8:
		return freq_hap<double>(xpMat, index_, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
NumericVector freq(XPtr<BigMatrix> pMat, const IntegerVector index_, const int threads = 0){

    if (threads == 0) {
        omp_set_num_threads(omp_get_num_procs());
    }else if(threads > 0) {
        omp_set_num_threads(threads);
    }

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();
	int j, k, m = index_.size();
	double p1 = 0.0, maf = 0.0;
	NumericVector freq(m);

	#pragma omp parallel for schedule(dynamic) private(j, p1, k, maf)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			p1 += bigm[index_[j]][k];
		}
		maf = p1 / (ind * 2);
		if(maf > 0.5){
			maf = 1 - maf;
		}
		freq[j] = maf;
	}
	return freq;
}

// [[Rcpp::export]]
NumericVector freq(SEXP pBigMat, const IntegerVector index_, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return freq<char>(xpMat, index_, threads);
	case 2:
		return freq<short>(xpMat, index_, threads);
	case 4:
		return freq<int>(xpMat, index_, threads);
	case 8:
		return freq<double>(xpMat, index_, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP BigStat(XPtr<BigMatrix> pMat, const IntegerVector index_, const int threads = 0){

    if (threads == 0) {
        omp_set_num_threads(omp_get_num_procs());
    }else if(threads > 0) {
        omp_set_num_threads(threads);
    }

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();
	int j, k, m = index_.size();
	double scale_mean, p1 = 0.0;
	NumericVector mean(m);
	NumericVector sd(m);
	NumericVector sum(m);

	#pragma omp parallel for private(p1, k)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			p1 += bigm[index_[j]][k];
		}
		sum[j] = p1;
		mean[j] = p1 / ind;
	}

	#pragma omp parallel for private(p1, k, scale_mean)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			scale_mean = (bigm[index_[j]][k] - mean[j]);
			p1 += scale_mean * scale_mean;
		}
		sd[j] = sqrt(p1);
	}
	return List::create(Named("mean") = mean, Named("sd") = sd, Named("sum") = sum);
}

// [[Rcpp::export]]
SEXP BigStat(SEXP pBigMat, const IntegerVector index_, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return BigStat<char>(xpMat, index_, threads);
	case 2:
		return BigStat<short>(xpMat, index_, threads);
	case 4:
		return BigStat<int>(xpMat, index_, threads);
	case 8:
		return BigStat<double>(xpMat, index_, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

// [[Rcpp::export]]
IntegerVector which_c(const NumericVector x, const double value, const int c = 1){

	IntegerVector eff_index_(x.size());
	int type_len = 0;
	bool logi;
	for(int i = 0; i < x.size(); i++){
		if(c == 1){
			logi = x[i] > value;
		}else if(c == 2){
			logi = x[i] >= value;
		}else if(c == 3){
			logi = x[i] < value;
		}else if(c == 4){
			logi = x[i] <= value;
		}else if(c == 5){
			logi = x[i] == value;
		}else if(c == 6){
			logi = (x[i] >= value) && (x[i] <= 1 - value);
		}else if(c == 7){
			logi = (x[i] < value) || (x[i] > 1 - value);
		}
		if(logi){
			eff_index_[type_len] = i;
			type_len++;
		}
	}
	IntegerVector eff_index(type_len);
	for(int i = 0; i < type_len; i++){
		eff_index[i] = eff_index_[i];
	}
	return eff_index;
}
