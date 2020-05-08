#include "stats.h"

template <typename T>
arma::uvec LDprune_c(XPtr<BigMatrix> pMat, const IntegerVector index, const double r2_cutoff = 0.2, const int threads=0, const bool verbose=true){
	
	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int m = index.size();
	int ind = pMat->nrow();
	int mm, i, j, k, findxi;
	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0;
	
	// MinimalProgressBar pb;

	IntegerVector index_ = index - 1;
	arma::vec freq_all = as<arma::vec>(freq(pMat, index_, threads));
	IntegerVector freq_sort_indx = as<Rcpp::IntegerVector>(wrap(sort_index(freq_all, "descend")));
	index_ = index_[freq_sort_indx];
	List Stat = BigStat(pMat, index_, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sd_all = Stat[1];
	NumericVector sum_all = Stat[2];

	// Progress p(m, verbose, pb);

	arma::uvec indx1(m); indx1.fill(true);
	arma::uvec indx2(m); indx2.fill(false);
	arma::uvec findx;

	for (j = 0; j < m; j++){
		if(indx1[j]){
			indx1[j] = 0;
			p1 = sd_all[j];
			m1 = mean_all[j];
			s1 = sum_all[j];
			findx = find(indx1);
			mm = findx.n_elem;
			indx2[freq_sort_indx[j]] = true;

			#pragma omp parallel for schedule(dynamic) private(s2, i, k, findxi, p12, p2, m2, r)
			for(i = 0; i < mm; i++){
				p12 = 0;
				findxi = findx[i];
				p2 = sd_all[findxi];
				m2 = mean_all[findxi];
				s2 = sum_all[findxi];
				for(k = 0; k < ind; k++){
					p12 += (genomat[index_[findxi]][k]) * (genomat[index_[j]][k]);
				}
				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
				r = p12 / (p1 * p2);
				if((r * r) > r2_cutoff)	indx1[findxi] = 0;
			}
		}else{
			continue;
		}
	}

	// arma::mat ldmat(m, m);

	// #pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
	// for (j = 0; j < m; j++){
	// 	ldmat(j, j) = 1;
	// 	p1 = sd_all[j];
	// 	m1 = mean_all[j];
	// 	s1 = sum_all[j];
	// 	for(i = j + 1; i < m; i++){
	// 		p12 = 0;
	// 		p2 = sd_all[i];
	// 		m2 = mean_all[i];
	// 		s2 = sum_all[i];
	// 		for(k = 0; k < ind; k++){
	// 			p12 += (genomat[index_[i]][k]) * (genomat[index_[j]][k]);
	// 		}
	// 		p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
	// 		r = p12 / (p1 * p2);
	// 		ldmat(j, i) = ldmat(i, j) = r * r;
	// 	}
	// }


	// for(uword i = 0; i < m; i++){
	// 	if(indx1[i]){
	// 		findx = find((ldmat.col(i) > r2_cutoff) && indx1);
	// 		indx1(findx).fill(false); 
	// 		indx2[freq_sort_indx[i]] = true;
	// 	}else{
	// 		continue;
	// 	}
	// }
	return indx2;
}

// [[Rcpp::export]]
arma::uvec LDprune_c(SEXP pBigMat, const IntegerVector index, const double r2_cutoff = 0.2, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return LDprune_c<char>(xpMat, index, r2_cutoff, threads, verbose);
	case 2:
		return LDprune_c<short>(xpMat, index, r2_cutoff, threads, verbose);
	case 4:
		return LDprune_c<int>(xpMat, index, r2_cutoff, threads, verbose);
	case 8:
		return LDprune_c<double>(xpMat, index, r2_cutoff, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
arma::uvec LDclump_c(XPtr<BigMatrix> pMat, const IntegerVector index, const double r2_cutoff, const arma::vec p, const int threads=0, const bool verbose=true){
	
	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int m = index.size();
	int ind = pMat->nrow();
	int mm, i, j, k, findxi;
	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0;
	
	// MinimalProgressBar pb;

	IntegerVector index_ = index - 1;
	IntegerVector freq_sort_indx = as<Rcpp::IntegerVector>(wrap(sort_index(p, "ascend")));
	index_ = index_[freq_sort_indx];
	List Stat = BigStat(pMat, index_, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sd_all = Stat[1];
	NumericVector sum_all = Stat[2];

	// Progress p(m, verbose, pb);

	arma::uvec indx1(m); indx1.fill(true);
	arma::uvec indx2(m); indx2.fill(false);
	arma::uvec findx;

	for (j = 0; j < m; j++){
		if(indx1[j]){
			indx1[j] = 0;
			p1 = sd_all[j];
			m1 = mean_all[j];
			s1 = sum_all[j];
			findx = find(indx1);
			mm = findx.n_elem;
			indx2[freq_sort_indx[j]] = true;

			#pragma omp parallel for schedule(dynamic) private(s2, i, k, findxi, p12, p2, m2, r)
			for(i = 0; i < mm; i++){
				p12 = 0;
				findxi = findx[i];
				p2 = sd_all[findxi];
				m2 = mean_all[findxi];
				s2 = sum_all[findxi];
				for(k = 0; k < ind; k++){
					p12 += (genomat[index_[findxi]][k]) * (genomat[index_[j]][k]);
				}
				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
				r = p12 / (p1 * p2);
				if((r * r) > r2_cutoff)	indx1[findxi] = 0;
			}
		}else{
			continue;
		}
	}
	return indx2;
}

// [[Rcpp::export]]
arma::uvec LDclump_c(SEXP pBigMat, const IntegerVector index, const double r2_cutoff, const arma::vec p, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return LDclump_c<char>(xpMat, index, r2_cutoff, p, threads, verbose);
	case 2:
		return LDclump_c<short>(xpMat, index, r2_cutoff, p, threads, verbose);
	case 4:
		return LDclump_c<int>(xpMat, index, r2_cutoff, p, threads, verbose);
	case 8:
		return LDclump_c<double>(xpMat, index, r2_cutoff, p, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
