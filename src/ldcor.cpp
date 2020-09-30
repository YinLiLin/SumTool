#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#include "stats.h"
#include "probar.h"
#include <string.h>

arma::rowvec cortest(arma::vec a, arma::vec b){
	arma::rowvec res(2);
	int n = a.n_elem;
	double r = as_scalar(cor(a, b));
	res[0] = r;
	double t = r * sqrt(n - 2) / sqrt(1 - r * r);
	double p = 2 * R::pt(abs(t), n - 2, false, false);
	res[1] = p;
	return res;
}

template <typename T>
arma::mat LDcor_c(XPtr<BigMatrix> pMat1, XPtr<BigMatrix> pMat2, const IntegerVector index1, const IntegerVector index2, const int threads=0){

	omp_setup(threads);

	MatrixAccessor<T> genomat1 = MatrixAccessor<T>(*pMat1);
	MatrixAccessor<T> genomat2 = MatrixAccessor<T>(*pMat2);

	int m = index1.size();
	int ind1 = pMat1->nrow();
	int ind2 = pMat2->nrow();
	int i, j, k;
	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0;
	List Stat1 = BigStat(pMat1, index1, threads);
	List Stat2 = BigStat(pMat2, index2, threads);
	NumericVector mean_all1 = Stat1[0];
	NumericVector sd_all1 = Stat1[1];
	NumericVector sum_all1 = Stat1[2];
	NumericVector mean_all2 = Stat2[0];
	NumericVector sd_all2 = Stat2[1];
	NumericVector sum_all2 = Stat2[2];

	arma::mat ldm1(m, m);
	arma::mat ldm2(m, m);
	arma::mat res(m, 2);

	#pragma omp parallel for private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
	for (j = 0; j < m; j++){
		p1 = sd_all1[j];
		m1 = mean_all1[j];
		s1 = sum_all1[j];
		for(i = j; i < m; i++){
			p12 = 0;
			p2 = sd_all1[i];
			m2 = mean_all1[i];
			s2 = sum_all1[i];
			for(k = 0; k < ind1; k++){
				p12 += (genomat1[index1[i]][k] * genomat1[index1[j]][k]);
			}
			p12 -= s1 * m2 + s2 * m1 - ind1 * m1 * m2;
			r = p12 / (p1 * p2);
			ldm1(i, j) = ldm1(j, i) = r;
		}
		p1 = sd_all2[j];
		m1 = mean_all2[j];
		s1 = sum_all2[j];
		for(i = j; i < m; i++){
			p12 = 0;
			p2 = sd_all2[i];
			m2 = mean_all2[i];
			s2 = sum_all2[i];
			for(k = 0; k < ind2; k++){
				p12 += (genomat2[index2[i]][k] * genomat2[index2[j]][k]);
			}
			p12 -= s1 * m2 + s2 * m1 - ind2 * m1 * m2;
			r = p12 / (p1 * p2);
			ldm2(i, j) = ldm2(j, i) = r;
		}
		res(j, 0) = as_scalar(cor(ldm1.col(j), ldm2.col(j)));
		double t = res(j, 0) * sqrt(m - 2) / sqrt(1 - res(j, 0) * res(j, 0));
		res(j, 1) = 2 * R::pt(abs(t), m - 2, false, false);
	}

	return res;
}

// [[Rcpp::export]]
arma::mat LDcor_c(SEXP pBigMat1, SEXP pBigMat2, const IntegerVector index1, const IntegerVector index2, const int threads=0){

	XPtr<BigMatrix> xpMat1(pBigMat1);
	XPtr<BigMatrix> xpMat2(pBigMat2);

	switch(xpMat1->matrix_type()){
	case 1:
		return LDcor_c<char>(xpMat1, xpMat2, index1, index2, threads);
	case 2:
		return LDcor_c<short>(xpMat1, xpMat2, index1, index2, threads);
	case 4:
		return LDcor_c<int>(xpMat1, xpMat2, index1, index2, threads);
	case 8:
		return LDcor_c<double>(xpMat1, xpMat2, index1, index2, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
