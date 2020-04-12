#include "stats.h"

template <typename T>
SEXP LDscore_c(XPtr<BigMatrix> pMat, const IntegerVector index, const bool r2 = true, const bool adjust = false, const int threads=0, const bool verbose=true){
	
	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int m = index.size();
	int ind = pMat->nrow();
	int i, j, k;
	double p1 = 0.0, q1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, q2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0, rr = 0.0;
	
	// MinimalProgressBar pb;

	IntegerVector index_ = index - 1;
	NumericVector freq_all = freq(pMat, index_, threads);
	List Stat = BigStat(pMat, index_, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sd_all = Stat[1];
	NumericVector sum_all = Stat[2];

	// Progress p(m, verbose, pb);

	arma::mat ldmat(m, m);

	if(r2){
		#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r, rr)
		for (j = 0; j < m; j++){
			ldmat(j, j) = 1;
			p1 = sd_all[j];
			m1 = mean_all[j];
			s1 = sum_all[j];
			for(i = j + 1; i < m; i++){
				p12 = 0;
				p2 = sd_all[i];
				m2 = mean_all[i];
				s2 = sum_all[i];
				for(k = 0; k < ind; k++){
					p12 += (genomat[index_[i]][k]) * (genomat[index_[j]][k]);
				}
				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
				// cout << j << "-" << p1 << "-" << m1 << "-" << p2 << "-" << m2 <<"-" << p12 << endl;
				r = p12 / (p1 * p2);
				if(adjust){
					rr = r * r;
					ldmat(j, i) = ldmat(i, j) = rr - ((1 - rr) / (ind - 2));
				}else{
					ldmat(j, i) = ldmat(i, j) = r * r;
				}
			}
		}
	}else{
		#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
		for (j = 0; j < m; j++){
			ldmat(j, j) = 1;
			p1 = sd_all[j];
			m1 = mean_all[j];
			s1 = sum_all[j];
			for(i = j + 1; i < m; i++){
				p12 = 0;
				p2 = sd_all[i];
				m2 = mean_all[i];
				s2 = sum_all[i];
				for(k = 0; k < ind; k++){
					p12 += (genomat[index_[i]][k]) * (genomat[index_[j]][k]);
				}
				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
				// cout << j << "-" << p1 << "-" << m1 << "-" << p2 << "-" << m2 <<"-" << p12 << endl;
				r = p12 / (p1 * p2);
				if(adjust){
					ldmat(j, i) = ldmat(i, j) = r - ((1 - r) / (ind - 2));
				}else{
					ldmat(j, i) = ldmat(i, j) = r;
				}
			}
		}
	}

	arma::mat ldscore(m, 4);
	for(int i = 0; i < m; i++){
		ldscore(i, 0) = freq_all[i];
		ldscore(i, 1) = m;
		ldscore(i, 3) = ldscore(i, 2) * m;
		ldscore(i, 2) = ldmat.col(i).mean();
	}

	return wrap(ldscore);
}

// [[Rcpp::export]]
SEXP LDscore_c(SEXP pBigMat, const IntegerVector index, const bool r2 = true, const bool adjust = false, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return LDscore_c<char>(xpMat, index, r2, adjust, threads, verbose);
	case 2:
		return LDscore_c<short>(xpMat, index, r2, adjust, threads, verbose);
	case 4:
		return LDscore_c<int>(xpMat, index, r2, adjust, threads, verbose);
	case 8:
		return LDscore_c<double>(xpMat, index, r2, adjust, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
