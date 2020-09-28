#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#include <RcppArmadillo.h>
#include "omp_set.h"
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

	omp_setup(threads);

    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    bool HasNA = false;
    int nc = pMat->ncol();
    int nr = pMat->nrow();

    #pragma omp parallel for schedule(dynamic) shared(HasNA)
    for (int j = 0; j < nc; j++) {
    	if(HasNA)	continue;
        for (int i = 0; i < nr; i++) {
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
NumericVector freq_snp(XPtr<BigMatrix> pMat, const IntegerVector index_, const int threads = 0){

    omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();
	int j, k, m = index_.size();
	int g;
	double p1 = 0.0;
	NumericVector freq(m);

	#pragma omp parallel for schedule(dynamic) private(j, p1, k, g)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			g = bigm[index_[j]][k];
			// p1 += (((g / 10) == 1) + ((g % 10) == 1));
			p1 += ((g == 11) ? (2) : (g == 22 ? (0) : (1)));
		}
		freq[j] = p1 / (ind * 2);
	}
	return freq;
}

// [[Rcpp::export]]
NumericVector freq_snp(SEXP pBigMat, const IntegerVector index_, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return freq_snp<char>(xpMat, index_, threads);
	case 2:
		return freq_snp<short>(xpMat, index_, threads);
	case 4:
		return freq_snp<int>(xpMat, index_, threads);
	case 8:
		return freq_snp<double>(xpMat, index_, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
double freq_hap(XPtr<BigMatrix> pMat, const int indx_1, const int indx_2){

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();
	int k;
	double p12 = 0.0;
	double freq;
	IntegerVector index_ = {indx_1, indx_2};

	NumericVector p = freq_snp(pMat, index_);
	for(k = 0; k < ind; k++){
		std::string g1 = std::to_string(bigm[indx_1][k]);
		std::string g2 = std::to_string(bigm[indx_2][k]);
		if(g1[0] == '1' && g2[0] == '1')	p12++;
		if(g1[1] == '1' && g2[1] == '1')	p12++;
	}
	p12 /= (ind * 2);
	freq = (p12 - p[0] * p[1]) / sqrt(p[0] * (1.0 - p[0]) * p[1] * (1.0 - p[1]));
	return freq;
}

// [[Rcpp::export]]
double freq_hap(SEXP pBigMat, const int indx_1, const int indx_2){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return freq_hap<char>(xpMat, indx_1, indx_2);
	case 2:
		return freq_hap<short>(xpMat, indx_1, indx_2);
	case 4:
		return freq_hap<int>(xpMat, indx_1, indx_2);
	case 8:
		return freq_hap<double>(xpMat, indx_1, indx_2);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
NumericVector freq_s(XPtr<BigMatrix> pMat, const IntegerVector index_, const int threads = 0){

    omp_setup(threads);

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
		freq[j] = maf;
	}
	return freq;
}

// [[Rcpp::export]]
NumericVector freq_s(SEXP pBigMat, const IntegerVector index_, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return freq_s<char>(xpMat, index_, threads);
	case 2:
		return freq_s<short>(xpMat, index_, threads);
	case 4:
		return freq_s<int>(xpMat, index_, threads);
	case 8:
		return freq_s<double>(xpMat, index_, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
double freq_h(XPtr<BigMatrix> pMat, const int indx_1, const int indx_2){

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();
	int k;
	double freq;

	double x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
	for(k = 0; k < ind; k++){
		int gi = bigm[indx_1][k];
		int gj = bigm[indx_2][k];
		if(gi == 0 && gj == 0) {
			x11 += 2.0;
		}
		else if (gi == 1 && gj == 1) {
			x11 += 1.0;
			x22 += 1.0;
		}
		else if (gi == 2 && gj == 2) {
			x22 += 2.0;
		}
		else if (gi == 0 && gj == 1) {
			x11 += 1.0;
			x12 += 1.0;
		}
		else if (gi == 1 && gj == 0) {
			x11 += 1.0;
			x21 += 1.0;
		}
		else if (gi == 0 && gj == 2) {
			x12 += 2.0;
		}
		else if (gi == 2 && gj == 0) {
			x21 += 2.0;
		}
		else if (gi == 1 && gj == 2) {
			x12 += 1.0;
			x22 += 1.0;
		} else {
			x21 += 1.0;
			x22 += 1.0;
		}
	}
	x11 /= ((double)(ind) * 2.0);
	x12 /= ((double)(ind) * 2.0);
	x21 /= ((double)(ind) * 2.0);
	x22 /= ((double)(ind) * 2.0);
	double p1 = x11 + x12;
	double p2 = x21 + x22;
	double q1 = x11 + x21;
	double q2 = x12 + x22;
	freq = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);

	return freq;
}

// [[Rcpp::export]]
double freq_h(SEXP pBigMat, const int indx_1, const int indx_2){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return freq_h<char>(xpMat, indx_1, indx_2);
	case 2:
		return freq_h<short>(xpMat, indx_1, indx_2);
	case 4:
		return freq_h<int>(xpMat, indx_1, indx_2);
	case 8:
		return freq_h<double>(xpMat, indx_1, indx_2);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP BigStat(XPtr<BigMatrix> pMat, const IntegerVector index_, const int threads = 0){

    omp_setup(threads);

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
	for(size_t i = 0; i < x.size(); i++){
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
