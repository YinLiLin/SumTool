#include "stats.h"

template <typename T>
SEXP sblup_bin(XPtr<BigMatrix> pMat, const double n_gwas, const IntegerVector typed_index, const arma::vec typed_value, const double lambda, const bool verbose = true, const int threads = 0){

	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}
	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int n = typed_index.size();
	int ind = pMat->nrow();
	int i, j, k;
	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, Dii = 0.0;
	IntegerVector index = typed_index - 1;

	List Stat = BigStat(pMat, index, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sd_all = Stat[1];
	NumericVector sum_all = Stat[2];

	arma::mat Rtt(n, n);
	arma::vec Dy(n);

	#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, Dii)
	for (j = 0; j < n; j++){
		p1 = sd_all[j];
		m1 = mean_all[j];
		s1 = sum_all[j];
		for(i = j + 1; i < n; i++){
			p12 = 0;
			p2 = sd_all[i];
			m2 = mean_all[i];
			s2 = sum_all[i];
			for(k = 0; k < ind; k++){
				p12 += (bigm[index[i]][k]) * (bigm[index[j]][k]);
			}
			p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
			// r = p12 / (p1 * p2);
			Rtt(i, j) = Rtt(j, i) = p12 * n_gwas / ind;
		}
		Dii = p1 * p1 * n_gwas / ind;
		Rtt(j, j) = Dii + lambda;
		Dy[j] = typed_value[j] * Dii;
	}

	arma::mat joint_eff = inv_sympd(Rtt) * Dy;
	// arma::vec eigval;
	// arma::mat eigvec;
	// eig_sym(eigval, eigvec, Rtt);
	
	// arma::mat Dy(typed_value.n_rows, typed_value.n_cols);
	// for(int i = 0; i < typed_value.n_rows; i++){
	// 	Dy.row(i) = typed_value.row(i) * Rtt(i, i);
	// }

	// arma::mat joint_eff(typed_value.n_rows, typed_value.n_cols);
	// for(i = 0; i < typed_value.n_cols; i++){
	// 	arma::rowvec eigveclambda = (1 / (eigval + lambda[i] / n_gwas[i])).t();
	// 	arma::mat eiglam = eigvec.each_row() % eigveclambda;
	// 	arma::mat iRtt = eiglam * eigvec.t();
	// 	joint_eff.col(i) = iRtt * Dy.col(i);
	// }

	return wrap(joint_eff);
}

// [[Rcpp::export]]
SEXP sblup_bin(SEXP pBigMat, const double n_gwas, const IntegerVector typed_index, const arma::vec typed_value, const double lambda, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return sblup_bin<char>(xpMat, n_gwas, typed_index, typed_value, lambda, verbose, threads);
	case 2:
		return sblup_bin<short>(xpMat, n_gwas, typed_index, typed_value, lambda, verbose, threads);
	case 4:
		return sblup_bin<int>(xpMat, n_gwas, typed_index, typed_value, lambda, verbose, threads);
	case 8:
		return sblup_bin<double>(xpMat, n_gwas, typed_index, typed_value, lambda, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
