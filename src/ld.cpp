#include "stats.h"
#include "probar.h"
#include <string.h>

template <typename T>
void SImpute_LD_bigm_c(arma::mat & genomat, XPtr<BigMatrix> ldMat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

	MatrixAccessor<T> ldmat = MatrixAccessor<T>(*ldMat);

	int m = genomat.n_cols;
	int ind = genomat.n_rows;
	int i, j, k, gi, gj;
	double p1 = 0.0, q1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, q2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0, x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
	bool sparse = false;
	MinimalProgressBar pb;
	IntegerVector index_all = seq(0, m - 1), index_;
	// NumericVector freq_all = freq(pMat, index_all, threads);
	// List Stat = BigStat(pMat, index_all, threads);
	NumericVector mean_all(m), sd_all(m), sum_all(m);
	if(!haps){
		#pragma omp parallel for schedule(dynamic)
		for(int i = 0; i < m; i++){
			mean_all[i] = mean(genomat.col(i));
			sum_all[i] = mean_all[i] * ind;
			sd_all[i] = sqrt(ind - 1) * stddev(genomat.col(i));
		}
	}
	if(chisq > 0)	sparse = true;

	Progress p(m, verbose, pb);

	if(index.isNotNull()){
		index_ = as<IntegerVector>(index) - 1;
		if(haps){

			#pragma omp parallel for schedule(dynamic) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
			for (j = 0; j < m; j++){
				if ( ! Progress::check_abort() ) {
					p.increment();
					for(i = 0; i < index_.size(); i++){
						if(index_[i] == j){
								r = 1 + lambda;
						}else{
							x11 = 0; x12 = 0; x21 = 0; x22 = 0;
							for(k = 0; k < ind; k++){
								gi = genomat(k, index_[i]);
								gj = genomat(k, index_[j]);
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
								}
								else if (gi == 2 && gj == 1) {
									x21 += 1.0;
									x22 += 1.0;
								}
							}
							x11 /= ((double)(ind) * 2.0);
							x12 /= ((double)(ind) * 2.0);
							x21 /= ((double)(ind) * 2.0);
							x22 /= ((double)(ind) * 2.0);
							p1 = x11 + x12;
							p2 = x21 + x22;
							q1 = x11 + x21;
							q2 = x12 + x22;
							r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
							if(sparse){
								if(r * r * ind < chisq)	r = 0.0;
							}
						}
						ldmat[i][j] = r;
					}
				}
			}
		}else{

			#pragma omp parallel for private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
			for (j = 0; j < m; j++){
				p1 = sd_all[j];
				m1 = mean_all[j];
				s1 = sum_all[j];
				if ( ! Progress::check_abort() ) {
					p.increment();

					for(i = 0; i < index_.size(); i++){
						if(index_[i] == j){
							r = 1 + lambda;
						}else{
							p12 = 0;
							p2 = sd_all[index_[i]];
							m2 = mean_all[index_[i]];
							s2 = sum_all[index_[i]];
							for(k = 0; k < ind; k++){
								p12 += (genomat(k, index_[i]) * genomat(k, j));
							}
							p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
							// cout << j << "-" << p1 << "-" << p2 << "-" << p12 << endl;
							r = p12 / (p1 * p2);
							if(sparse){
								if(r * r * ind < chisq)	r = 0.0;
							}
						}
						ldmat[i][j] = r;
					}
				}
			}
		}
	}else{
		if(haps){

			#pragma omp parallel for schedule(dynamic) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
			for (j = 0; j < m; j++){
				ldmat[j][j] = 1 + lambda;
				if ( ! Progress::check_abort() ) {
					p.increment();
					for(i = j + 1; i < m; i++){
						x11 = 0; x12 = 0; x21 = 0; x22 = 0;
						for(k = 0; k < ind; k++){
							gi = genomat(k, i);
							gj = genomat(k, j);
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
							}
							else if (gi == 2 && gj == 1) {
								x21 += 1.0;
								x22 += 1.0;
							}
						}
						x11 /= ((double)(ind) * 2.0);
						x12 /= ((double)(ind) * 2.0);
						x21 /= ((double)(ind) * 2.0);
						x22 /= ((double)(ind) * 2.0);
						p1 = x11 + x12;
						p2 = x21 + x22;
						q1 = x11 + x21;
						q2 = x12 + x22;
						r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
						if(sparse){
							if(r * r * ind < chisq)	r = 0.0;
						}
						ldmat[i][j] = ldmat[j][i] = r;
					}
				}
			}
		}else{

			#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
			for (j = 0; j < m; j++){
				ldmat[j][j] = 1 + lambda;
				p1 = sd_all[j];
				m1 = mean_all[j];
				s1 = sum_all[j];
				if ( ! Progress::check_abort() ) {
					p.increment();

					for(i = j + 1; i < m; i++){
						p12 = 0;
						p2 = sd_all[i];
						m2 = mean_all[i];
						s2 = sum_all[i];
						for(k = 0; k < ind; k++){
							p12 += (genomat(k, i) * genomat(k, j));
						}
						p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
						// cout << j << "-" << p1 << "-" << m1 << "-" << p2 << "-" << m2 <<"-" << p12 << endl;
						r = p12 / (p1 * p2);
						if(sparse){
							if(r * r * ind < chisq)	r = 0.0;
						}
						ldmat[i][j] = ldmat[j][i] = r;
					}
				}
			}
		}
	}
	return;
}

// [[Rcpp::export]]
void SImpute_LD_bigm_c(arma::mat & genomat, SEXP ldbigmat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(ldbigmat);

	switch(xpMat->matrix_type()){
	case 1:
		return SImpute_LD_bigm_c<char>(genomat, xpMat, index, chisq, lambda, haps, threads, verbose);
	case 2:
		return SImpute_LD_bigm_c<short>(genomat, xpMat, index, chisq, lambda, haps, threads, verbose);
	case 4:
		return SImpute_LD_bigm_c<int>(genomat, xpMat, index, chisq, lambda, haps, threads, verbose);
	case 8:
		return SImpute_LD_bigm_c<double>(genomat, xpMat, index, chisq, lambda, haps, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

// template <typename T>
// void SImpute_LD_bigm_c(XPtr<BigMatrix> pMat, const SEXP ldbigmat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

// 	if (threads == 0) {
// 		omp_set_num_threads(omp_get_num_procs());
// 	}else if(threads > 0) {
// 		omp_set_num_threads(threads);
// 	}

// 	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);
// 	XPtr<BigMatrix> ldbigMat(ldbigmat);
// 	MatrixAccessor<double> ldmat(*ldbigMat);

// 	int m = pMat->ncol();
// 	int ind = pMat->nrow();
// 	int i, j, k, gi, gj;
// 	double p1 = 0.0, q1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, q2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0, x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
// 	bool sparse = false;
// 	MinimalProgressBar pb;
// 	IntegerVector index_all = seq(0, m - 1), index_;
// 	// NumericVector freq_all = freq(pMat, index_all, threads);
// 	List Stat = BigStat(pMat, index_all, threads);
// 	NumericVector mean_all, sd_all, sum_all;
// 	if(!haps){
// 		mean_all = Stat[0];
// 		sd_all = Stat[1];
// 		sum_all = Stat[2];
// 	}
// 	if(chisq > 0)	sparse = true;

// 	Progress p(m, verbose, pb);

// 	if(index.isNotNull()){
// 		index_ = as<IntegerVector>(index) - 1;
// 		arma::vec Rit(index_.size());
// 		if(haps){
// 			for (j = 0; j < m; j++){
// 				if ( ! Progress::check_abort() ) {
// 					p.increment();

// 					#pragma omp parallel for schedule(dynamic) private(i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
// 					for(i = 0; i < index_.size(); i++){
// 						if(index_[i] == j){
// 								r = 1 + lambda;
// 						}else{
// 							x11 = 0; x12 = 0; x21 = 0; x22 = 0;
// 							for(k = 0; k < ind; k++){
// 								gi = genomat[index_[i]][k];
// 								gj = genomat[index_[j]][k];
// 								if(gi == 0 && gj == 0) {
// 									x11 += 2.0;
// 								}
// 								else if (gi == 1 && gj == 1) {
// 									x11 += 1.0;
// 									x22 += 1.0;
// 								}
// 								else if (gi == 2 && gj == 2) {
// 									x22 += 2.0;
// 								}
// 								else if (gi == 0 && gj == 1) {
// 									x11 += 1.0;
// 									x12 += 1.0;
// 								}
// 								else if (gi == 1 && gj == 0) {
// 									x11 += 1.0;
// 									x21 += 1.0;
// 								}
// 								else if (gi == 0 && gj == 2) {
// 									x12 += 2.0;
// 								}
// 								else if (gi == 2 && gj == 0) {
// 									x21 += 2.0;
// 								}
// 								else if (gi == 1 && gj == 2) {
// 									x12 += 1.0;
// 									x22 += 1.0;
// 								}
// 								else if (gi == 2 && gj == 1) {
// 									x21 += 1.0;
// 									x22 += 1.0;
// 								}
// 							}
// 							x11 /= ((double)(ind) * 2.0);
// 							x12 /= ((double)(ind) * 2.0);
// 							x21 /= ((double)(ind) * 2.0);
// 							x22 /= ((double)(ind) * 2.0);
// 							p1 = x11 + x12;
// 							p2 = x21 + x22;
// 							q1 = x11 + x21;
// 							q2 = x12 + x22;
// 							r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
// 							if(sparse){
// 								if(r * r * ind < chisq)	r = 0.0;
// 							}
// 						}
// 						Rit[i] = r;
// 					}
// 					for(int l = 0; l < index_.size(); l++){
// 						ldmat[l][j] = Rit[l];
// 					}
// 				}
// 			}
// 		}else{
// 			for (j = 0; j < m; j++){
// 				p1 = sd_all[j];
// 				m1 = mean_all[j];
// 				s1 = sum_all[j];
// 				if ( ! Progress::check_abort() ) {
// 					p.increment();

// 					#pragma omp parallel for private(i, s2, k, p12, p2, m2, r)
// 					for(i = 0; i < index_.size(); i++){
// 						if(index_[i] == j){
// 							r = 1 + lambda;
// 						}else{
// 							p12 = 0;
// 							p2 = sd_all[index_[i]];
// 							m2 = mean_all[index_[i]];
// 							s2 = sum_all[index_[i]];
// 							for(k = 0; k < ind; k++){
// 								p12 += (genomat[index_[i]][k]) * (genomat[j][k]);
// 							}
// 							p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
// 							// cout << j << "-" << p1 << "-" << p2 << "-" << p12 << endl;
// 							r = p12 / (p1 * p2);
// 							if(sparse){
// 								if(r * r * ind < chisq)	r = 0.0;
// 							}
// 						}
// 						Rit[i] = r;
// 					}
// 					for(int l = 0; l < index_.size(); l++){
// 						ldmat[l][j] = Rit[l];
// 					}
// 				}
// 			}
// 		}
// 	}else{
// 		arma::vec Rit(m);
// 		if(haps){
// 			for (j = 0; j < m; j++){
// 				ldmat[j][j] = 1 + lambda;
// 				if ( ! Progress::check_abort() ) {
// 					p.increment();

// 					#pragma omp parallel for schedule(dynamic) private(i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
// 					for(i = j + 1; i < m; i++){
// 						x11 = 0; x12 = 0; x21 = 0; x22 = 0;
// 						for(k = 0; k < ind; k++){
// 							gi = genomat[i][k];
// 							gj = genomat[j][k];
// 							if(gi == 0 && gj == 0) {
// 								x11 += 2.0;
// 							}
// 							else if (gi == 1 && gj == 1) {
// 								x11 += 1.0;
// 								x22 += 1.0;
// 							}
// 							else if (gi == 2 && gj == 2) {
// 								x22 += 2.0;
// 							}
// 							else if (gi == 0 && gj == 1) {
// 								x11 += 1.0;
// 								x12 += 1.0;
// 							}
// 							else if (gi == 1 && gj == 0) {
// 								x11 += 1.0;
// 								x21 += 1.0;
// 							}
// 							else if (gi == 0 && gj == 2) {
// 								x12 += 2.0;
// 							}
// 							else if (gi == 2 && gj == 0) {
// 								x21 += 2.0;
// 							}
// 							else if (gi == 1 && gj == 2) {
// 								x12 += 1.0;
// 								x22 += 1.0;
// 							}
// 							else if (gi == 2 && gj == 1) {
// 								x21 += 1.0;
// 								x22 += 1.0;
// 							}
// 						}
// 						x11 /= ((double)(ind) * 2.0);
// 						x12 /= ((double)(ind) * 2.0);
// 						x21 /= ((double)(ind) * 2.0);
// 						x22 /= ((double)(ind) * 2.0);
// 						p1 = x11 + x12;
// 						p2 = x21 + x22;
// 						q1 = x11 + x21;
// 						q2 = x12 + x22;
// 						r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
// 						if(sparse){
// 							if(r * r * ind < chisq)	r = 0.0;
// 						}
// 						Rit[i] = r;
// 					}
// 					for(int l = j + 1; l < m; l++){
// 						ldmat[l][j] = ldmat[j][l] = Rit[l];
// 					}
// 				}
// 			}
// 		}else{
// 			for (j = 0; j < m; j++){
// 				ldmat[j][j] = 1 + lambda;
// 				p1 = sd_all[j];
// 				m1 = mean_all[j];
// 				s1 = sum_all[j];
// 				if ( ! Progress::check_abort() ) {
// 					p.increment();

// 					#pragma omp parallel for private(i, s2, k, p12, p2, m2, r)
// 					for(i = j + 1; i < m; i++){
// 						p12 = 0;
// 						p2 = sd_all[i];
// 						m2 = mean_all[i];
// 						s2 = sum_all[i];
// 						for(k = 0; k < ind; k++){
// 							p12 += (genomat[i][k]) * (genomat[j][k]);
// 						}
// 						p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
// 						// cout << j << "-" << p1 << "-" << p2 << "-" << p12 << endl;
// 						r = p12 / (p1 * p2);
// 						if(sparse){
// 							if(r * r * ind < chisq)	r = 0.0;
// 						}
// 						Rit[i] = r;
// 					}
// 					for(int l = j + 1; l < m; l++){
// 						ldmat[l][j] = ldmat[j][l] = Rit[l];
// 					}
// 				}
// 			}
// 		}
// 	}
// 	return;
// }

// // [[Rcpp::export]]
// void SImpute_LD_bigm_c(SEXP pBigMat, SEXP ldbigmat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

// 	XPtr<BigMatrix> xpMat(pBigMat);

// 	switch(xpMat->matrix_type()){
// 	case 1:
// 		return SImpute_LD_bigm_c<char>(xpMat, ldbigmat, index, chisq, lambda, haps, threads, verbose);
// 	case 2:
// 		return SImpute_LD_bigm_c<short>(xpMat, ldbigmat, index, chisq, lambda, haps, threads, verbose);
// 	case 4:
// 		return SImpute_LD_bigm_c<int>(xpMat, ldbigmat, index, chisq, lambda, haps, threads, verbose);
// 	case 8:
// 		return SImpute_LD_bigm_c<double>(xpMat, ldbigmat, index, chisq, lambda, haps, threads, verbose);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

template <typename T>
SEXP SImpute_LD_norm_c(XPtr<BigMatrix> pMat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int m = pMat->ncol(), mc;
	int ind = pMat->nrow();
	bool sparse = false;
	int i, j, k, gi, gj;
	double p1 = 0.0, q1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, q2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0, x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
	
	MinimalProgressBar pb;
	IntegerVector index_all = seq(0, m - 1), index_;
	// NumericVector freq_all = freq(pMat, index_all, threads);
	NumericVector mean_all, sd_all, sum_all;
	if(!haps){
		List Stat = BigStat(pMat, index_all, threads);
		mean_all = Stat[0];
		sd_all = Stat[1];
		sum_all = Stat[2];
	}

	if(chisq > 0)	sparse = true;

	Progress p(m, verbose, pb);

	if(index.isNotNull()){
		index_ = as<IntegerVector>(index) - 1;
		mc = index_.size();
	}else{
		mc = m;
	}

	arma::mat ldmat(m, mc);

	if(index.isNotNull()){
		if(haps){

			#pragma omp parallel for schedule(dynamic) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
			for (j = 0; j < m; j++){
				if ( ! Progress::check_abort() ) {
					p.increment();
					for(i = 0; i < mc; i++){
						if(index_[i] == j){
								r = 1 + lambda;
						}else{
							x11 = 0; x12 = 0; x21 = 0; x22 = 0;
							for(k = 0; k < ind; k++){
								gi = genomat[index_[i]][k];
								gj = genomat[index_[j]][k];
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
								}
								else if (gi == 2 && gj == 1) {
									x21 += 1.0;
									x22 += 1.0;
								}
							}
							x11 /= ((double)(ind) * 2.0);
							x12 /= ((double)(ind) * 2.0);
							x21 /= ((double)(ind) * 2.0);
							x22 /= ((double)(ind) * 2.0);
							p1 = x11 + x12;
							p2 = x21 + x22;
							q1 = x11 + x21;
							q2 = x12 + x22;
							r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
							if(sparse){
								if(r * r * ind < chisq)	r = 0.0;
							}
						}
						ldmat(j, i) = r;
					}
				}
			}
		}else{

			#pragma omp parallel for private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
			for (j = 0; j < m; j++){
				p1 = sd_all[j];
				m1 = mean_all[j];
				s1 = sum_all[j];
				if ( ! Progress::check_abort() ) {
					p.increment();

					for(i = 0; i < mc; i++){
						if(index_[i] == j){
							r = 1 + lambda;
						}else{
							p12 = 0;
							p2 = sd_all[index_[i]];
							m2 = mean_all[index_[i]];
							s2 = sum_all[index_[i]];
							for(k = 0; k < ind; k++){
								p12 += (genomat[index_[i]][k]) * (genomat[j][k]);
							}
							p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
							// cout << j << "-" << p1 << "-" << p2 << "-" << p12 << endl;
							r = p12 / (p1 * p2);
							if(sparse){
								if(r * r * ind < chisq)	r = 0.0;
							}
						}
						ldmat(j, i) = r;
					}
				}
			}
		}
	}else{
		if(haps){

			#pragma omp parallel for schedule(dynamic) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
			for (j = 0; j < m; j++){
				ldmat(j, j) = 1 + lambda;
				if ( ! Progress::check_abort() ) {
					p.increment();
					for(i = j + 1; i < m; i++){
						x11 = 0; x12 = 0; x21 = 0; x22 = 0;
						for(k = 0; k < ind; k++){
							gi = genomat[i][k];
							gj = genomat[j][k];
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
							}
							else if (gi == 2 && gj == 1) {
								x21 += 1.0;
								x22 += 1.0;
							}
						}
						x11 /= ((double)(ind) * 2.0);
						x12 /= ((double)(ind) * 2.0);
						x21 /= ((double)(ind) * 2.0);
						x22 /= ((double)(ind) * 2.0);
						p1 = x11 + x12;
						p2 = x21 + x22;
						q1 = x11 + x21;
						q2 = x12 + x22;
						r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
						if(sparse){
							if(r * r * ind < chisq)	r = 0.0;
						}
						ldmat(i, j) = ldmat(j, i) = r;
					}
				}
			}
		}else{

			#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
			for (j = 0; j < m; j++){
				ldmat(j, j) = 1 + lambda;
				p1 = sd_all[j];
				m1 = mean_all[j];
				s1 = sum_all[j];
				if ( ! Progress::check_abort() ) {
					p.increment();

					for(i = j + 1; i < m; i++){
						p12 = 0;
						p2 = sd_all[i];
						m2 = mean_all[i];
						s2 = sum_all[i];
						for(k = 0; k < ind; k++){
							p12 += (genomat[i][k]) * (genomat[j][k]);
						}
						p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
						// cout << j << "-" << p1 << "-" << m1 << "-" << p2 << "-" << m2 <<"-" << p12 << endl;
						r = p12 / (p1 * p2);
						if(sparse){
							if(r * r * ind < chisq)	r = 0.0;
						}
						ldmat(j, i) = ldmat(i, j) = r;
					}
				}
			}
		}
	}
	return wrap(ldmat);
}

// [[Rcpp::export]]
SEXP SImpute_LD_norm_c(SEXP pBigMat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return SImpute_LD_norm_c<char>(xpMat, index, chisq, lambda, haps, threads, verbose);
	case 2:
		return SImpute_LD_norm_c<short>(xpMat, index, chisq, lambda, haps, threads, verbose);
	case 4:
		return SImpute_LD_norm_c<int>(xpMat, index, chisq, lambda, haps, threads, verbose);
	case 8:
		return SImpute_LD_norm_c<double>(xpMat, index, chisq, lambda, haps, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP SImpute_LD_sparse_c(XPtr<BigMatrix> pMat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int m = pMat->ncol(), mc;
	int ind = pMat->nrow();
	bool sparse = false;
	int i, j, k, gi, gj;
	double p1 = 0.0, q1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, q2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0, x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
	
	MinimalProgressBar pb;
	IntegerVector index_all = seq(0, m - 1), index_;
	// NumericVector freq_all = freq(pMat, index_all, threads);
	NumericVector mean_all, sd_all, sum_all;
	if(!haps){
		List Stat = BigStat(pMat, index_all, threads);
		mean_all = Stat[0];
		sd_all = Stat[1];
		sum_all = Stat[2];
	}

	if(chisq > 0)	sparse = true;

	Progress p(m, verbose, pb);

	if(index.isNotNull()){
		index_ = as<IntegerVector>(index) - 1;
		mc = index_.size();
	}else{
		mc = m;
	}

	std::vector<std::string> item;
	List res(m * m);

	#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, s2, i, k, p12, p2, m2, r)
	for (j = 0; j < m; j++){
		// ldmat(j, j) = 1 + lambda;
		p1 = sd_all[j];
		m1 = mean_all[j];
		s1 = sum_all[j];
		if ( ! Progress::check_abort() ) {
			p.increment();

			for(i = j + 1; i < m; i++){
				p12 = 0;
				p2 = sd_all[i];
				m2 = mean_all[i];
				s2 = sum_all[i];
				for(k = 0; k < ind; k++){
					p12 += (genomat[i][k]) * (genomat[j][k]);
				}
				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
				// cout << j << "-" << p1 << "-" << m1 << "-" << p2 << "-" << m2 <<"-" << p12 << endl;
				r = p12 / (p1 * p2);
				if(sparse){
					if(r * r * ind < chisq)	r = 0.0;
				}
				if(r){
					std::string elem;
					elem = to_string(j) + " " + to_string(i) + " " + to_string(r);
					// item.push_back(elem);
					res[j * m + i] = elem;
				}
				// ldmat(j, i) = ldmat(i, j) = r;
			}
		}
	}
	return wrap(res);
}

// [[Rcpp::export]]
SEXP SImpute_LD_sparse_c(SEXP pBigMat, const Nullable<IntegerVector> index = R_NilValue, const int chisq = 0, const double lambda = 0, const bool haps = false, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return SImpute_LD_sparse_c<char>(xpMat, index, chisq, lambda, haps, threads, verbose);
	case 2:
		return SImpute_LD_sparse_c<short>(xpMat, index, chisq, lambda, haps, threads, verbose);
	case 4:
		return SImpute_LD_sparse_c<int>(xpMat, index, chisq, lambda, haps, threads, verbose);
	case 8:
		return SImpute_LD_sparse_c<double>(xpMat, index, chisq, lambda, haps, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
