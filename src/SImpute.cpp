#include "stats.h"
#include "probar.h"

template <typename T>
SEXP SImputeZ_bin_c(XPtr<BigMatrix> pMat, const IntegerVector typed_index, const arma::mat typed_value, const double lambda = 0.0, const double maf = 0.01, const bool haps = false, const bool verbose = true, const int threads = 0){

	omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int n = pMat->ncol();
	int ind = pMat->nrow();
	int m, i, j, k, gi, gj;

	double p1 = 0.0, q1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, q2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0, x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
	// arma::mat typed_value = Rcpp::as<arma::mat>(typed_value_R);
	IntegerVector index = seq(0, n - 1);
	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
	// arma::mat r2pred(n, typed_value.n_cols);
	MinimalProgressBar pb;
	// MinimalProgressBar_plus pb;
	IntegerVector typed_index_ = typed_index - 1;
	IntegerVector nontyped_index_ = rev(Rcpp::setdiff(index, typed_index_));
	NumericVector mean_all, sd_all, sum_all, typed_sd, typed_mean, typed_sum;

	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index_.size() << endl;
	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

	if(verbose)	Rcerr << "MAF calculating" << endl;
	NumericVector freq_all = freq(pMat, index, threads);
	if(!haps){
		List Stat = BigStat(pMat, index, threads);
		mean_all = Stat[0];
		sd_all = Stat[1];
		sum_all = Stat[2];
	}

	IntegerVector eff_index = which_c(freq_all[typed_index_], maf, 6);

	int typed_n = typed_index_.size();
	int typed_ncol = typed_value.n_cols;
	for(i = 0; i < typed_n; i++){
		for(j = 0; j < typed_ncol; j++){
			impz(typed_index_[i], j) = typed_value(i, j);
			// r2pred[typed_index_[i], j] = 1;
		}
	}

	IntegerVector index_ = typed_index_[eff_index];
	arma::mat eff_typed = typed_value.rows(as<uvec>(eff_index));
	m = index_.size();

	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;

	arma::mat Rtt(m, m);

	if(haps){

		#pragma omp parallel for schedule(dynamic) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
		for (j = 0; j < m; j++){
			for(i = j + 1; i < m; i++){
				x11 = 0; x12 = 0; x21 = 0; x22 = 0;
				for(k = 0; k < ind; k++){
					gi = bigm[index_[i]][k];
					gj = bigm[index_[j]][k];
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
				p1 = x11 + x12;
				p2 = x21 + x22;
				q1 = x11 + x21;
				q2 = x12 + x22;
				r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
				Rtt(i, j) = Rtt(j, i) = r;
			}
			Rtt(j, j) = 1 + lambda;
		}
	}else{
		typed_sd = sd_all[index_];
		typed_mean = mean_all[index_];
		typed_sum = sum_all[index_];

		#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r)
		for (j = 0; j < m; j++){
			p1 = typed_sd[j];
			m1 = typed_mean[j];
			s1 = typed_sum[j];
			for(i = j + 1; i < m; i++){
				p12 = 0;
				p2 = typed_sd[i];
				m2 = typed_mean[i];
				s2 = typed_sum[i];
				for(k = 0; k < ind; k++){
					p12 += (bigm[index_[i]][k]) * (bigm[index_[j]][k]);
				}
				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
				r = p12 / (p1 * p2);
				Rtt(i, j) = Rtt(j, i) = r;
			}
			Rtt(j, j) = 1 + lambda;
		}
	}

	if(verbose)	Rcerr << "Calculating the inverse of sigma_tt matrix" << endl;
	arma::mat iRtt = inv_sympd(Rtt);

	// arma::mat iRtt = inv(Rtt);

	// if(chisq > 0){
	// 	sparse = true;
	// 	if(verbose)	Rcerr << "Sparse sigma_it matrix generates at (Chisq): " << chisq << endl;
	// }

	// arma::vec Rit(m);
	if(verbose)	Rcerr << "Imputing started..." << endl;
	arma::mat beta;

	// Progress p(nontyped_index_.size(), verbose);
	int nontyped_n = nontyped_index_.size();
	Progress p(nontyped_n, verbose, pb);
	arma::vec Rit(m); Rit.fill(0);

	if(haps){

		#pragma omp parallel for schedule(dynamic) firstprivate(Rit) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r, beta)
		for(j = 0; j < nontyped_n; j++){
			p1 = freq_all[nontyped_index_[j]];
			if ( ! Progress::check_abort() ) {
				p.increment();
				if(p1 < maf){
					// impz[nontyped_index_[j]] = 0;
					// r2pred[nontyped_index_[j]] = 0;
				}else{
					for(i = 0; i < m; i++){
						x11 = 0; x12 = 0; x21 = 0; x22 = 0;
						for(k = 0; k < ind; k++){
							gi = bigm[index_[i]][k];
							gj = bigm[nontyped_index_[j]][k];
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
						p1 = x11 + x12;
						p2 = x21 + x22;
						q1 = x11 + x21;
						q2 = x12 + x22;
						r = (x11 * x22 - x12 * x21) / (sqrt(p1 * p2 * q1 * q2));
						Rit[i] = r;
					}
					beta = Rit.t() * iRtt;
					impz.row(nontyped_index_[j]) = beta * eff_typed;
				}
			}
		}
	}else{

		#pragma omp parallel for schedule(dynamic) firstprivate(Rit) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r, beta)
		for(j = 0; j < nontyped_n; j++){

			p1 = freq_all[nontyped_index_[j]];
			m1 = mean_all[nontyped_index_[j]];
			s1 = sum_all[nontyped_index_[j]];

			if ( ! Progress::check_abort() ) {
				p.increment();
				if(p1 < maf){
					// not impute
					// impz[nontyped_index_[j]] = 0;
					// r2pred[nontyped_index_[j]] = 0;
				}else{

					// impute
					p1 = sd_all[nontyped_index_[j]];
					for(i = 0; i < m; i++){
						p12 = 0;
						p2 = typed_sd[i];
						m2 = typed_mean[i];
						s2 = typed_sum[i];
						for(k = 0; k < ind; k++){
							p12 += (bigm[index_[i]][k]) * (bigm[nontyped_index_[j]][k]);
						}
						p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
						r = p12 / (p1 * p2);
						// if(sparse){
						// 	if(r * r * ind < chisq){
						// 		r = 0.0;
						// 	}
						// }
						Rit[i] = r;
					}

					// get beta
					beta = Rit.t() * iRtt;
					// impute missing and get var
					// impz[nontyped_index_[j], 1] = arma::as_scalar(beta * eff_typed);
					impz.row(nontyped_index_[j]) = beta * eff_typed;
					// r2pred.row(nontyped_index_[j]) = beta * Rit;
				}
			}
		}
	}
	
	// return List::create(Named("impz") = impz, Named("r2pred") = r2pred);
	return wrap(impz);
}

// [[Rcpp::export]]
SEXP SImputeZ_bin_c(SEXP pBigMat, const IntegerVector typed_index, const arma::mat typed_value, const double lambda = 0.0, const double maf = 0.01, const bool haps = false, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return SImputeZ_bin_c<char>(xpMat, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	case 2:
		return SImputeZ_bin_c<short>(xpMat, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	case 4:
		return SImputeZ_bin_c<int>(xpMat, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	case 8:
		return SImputeZ_bin_c<double>(xpMat, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP SImputeZ_ld_bin_c(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pMat_G, const IntegerVector typed_index, const arma::mat typed_value, const double lambda = 0.0, const double maf = 0.01, const bool haps = false, const bool verbose = true, const int threads = 0){

	omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
	MatrixAccessor<T> bigm_G = MatrixAccessor<T>(*pMat_G);

	int n = pMat->ncol();
	int n_gwas = pMat_G->ncol();
	int ind = pMat->nrow();
	int ind_gwas = pMat_G->nrow();
	int m, i, j, k, gi, gj;
	double p1 = 0.0, q1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, q2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0, x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
	// arma::mat typed_value = Rcpp::as<arma::mat>(typed_value_R);
	IntegerVector index = seq(0, n - 1);
	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
	// arma::mat r2pred(n, typed_value.n_cols);
	MinimalProgressBar pb;
	// MinimalProgressBar_plus pb;
	IntegerVector typed_index_ = typed_index - 1;
	IntegerVector nontyped_index_ = rev(Rcpp::setdiff(index, typed_index_));
	NumericVector mean_all, sd_all, sum_all, mean_gwas, sd_gwas, typed_sd, sum_gwas, typed_mean, typed_sum;

	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index_.size() << endl;
	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

	if(verbose)	Rcerr << "MAF calculating" << endl;
	NumericVector freq_all = freq(pMat, index, threads);
	IntegerVector eff_index = which_c(freq_all[typed_index_], maf, 6);
	if(!haps){
		List Stat = BigStat(pMat, index, threads);
		mean_all = Stat[0];
		sd_all = Stat[1];
		sum_all = Stat[2];
		List Stat_gwas = BigStat(pMat_G, seq(0, n_gwas - 1), threads);
		mean_gwas = Stat_gwas[0];
		sd_gwas = Stat_gwas[1];
		sum_gwas = Stat_gwas[2];
		typed_sd = sd_gwas[eff_index];
		typed_mean = mean_gwas[eff_index];
		typed_sum = sum_gwas[eff_index];
	}

	// typed SNPs needn't to impute
	int typed_n = typed_index_.size();
	int typed_ncol = typed_value.n_cols;
	for(i = 0; i < typed_n; i++){
		for(j = 0; j < typed_ncol; j++){
			impz(typed_index_[i], j) = typed_value(i, j);
			// r2pred[typed_index_[i], j] = 1;
		}
	}

	IntegerVector index_ = typed_index_[eff_index];
	int eff_index_n = eff_index.size();
	arma::mat eff_typed(eff_index_n, typed_value.n_cols);
	for(i = 0; i < eff_index_n; i++){
		eff_typed.row(i) = typed_value.row(eff_index[i]);
	}
	m = index_.size();

	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;

	arma::mat Rtt(m, m);

	if(haps){

		#pragma omp parallel for schedule(dynamic) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r)
		for (j = 0; j < m; j++){
			for(i = j + 1; i < m; i++){
				x11 = 0; x12 = 0; x21 = 0; x22 = 0;
				for(k = 0; k < ind_gwas; k++){
					gi = bigm_G[eff_index[i]][k];
					gj = bigm_G[eff_index[j]][k];
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
				x11 /= ((double)(ind_gwas) * 2.0);
				x12 /= ((double)(ind_gwas) * 2.0);
				x21 /= ((double)(ind_gwas) * 2.0);
				x22 /= ((double)(ind_gwas) * 2.0);
				p1 = x11 + x12;
				p2 = x21 + x22;
				q1 = x11 + x21;
				q2 = x12 + x22;
				r = (x11 * x22 - x12 * x21) / sqrt(p1 * p2 * q1 * q2);
				Rtt(i, j) = Rtt(j, i) = r;
			}
			Rtt(j, j) = 1 + lambda;
		}
	}else{

		#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r)
		for (j = 0; j < m; j++){
			p1 = typed_sd[j];
			m1 = typed_mean[j];
			s1 = typed_sum[j];
			for(i = j + 1; i < m; i++){
				p12 = 0;
				p2 = typed_sd[i];
				m2 = typed_mean[i];
				s2 = typed_sum[i];
				for(k = 0; k < ind_gwas; k++){
					p12 += (bigm_G[eff_index[i]][k]) * (bigm_G[eff_index[j]][k]);
				}
				p12 -= s1 * m2 + s2 * m1 - ind_gwas * m1 * m2;
				r = p12 / (p1 * p2);
				Rtt(i, j) = Rtt(j, i) = r;
			}
			Rtt(j, j) = 1 + lambda;
		}
	}

	// Rcout << Rtt(0, 0) << Rtt(1, 0) << Rtt(2, 0) << Rtt(3, 0) << Rtt(4, 0) << endl;
	if(verbose)	Rcerr << "Calculating the inverse of sigma_tt matrix" << endl;
	arma::mat iRtt = inv_sympd(Rtt);

	// arma::mat iRtt = inv(Rtt);

	// if(chisq > 0){
	// 	sparse = true;
	// 	if(verbose)	Rcerr << "Sparse sigma_it matrix generates at (Chisq): " << chisq << endl;
	// }

	// arma::vec Rit(m);
	if(verbose)	Rcerr << "Imputing started..." << endl;
	arma::mat beta;
	arma::vec Rit(m); Rit.fill(0);

	int nontyped_n = nontyped_index_.size();
	Progress p(nontyped_n, verbose, pb);

	if(haps){

		#pragma omp parallel for schedule(dynamic) firstprivate(Rit) private(j, i, x11, x12, x21, x22, k, gi, gj, p1, p2, q1, q2, r, beta)
		for(j = 0; j < nontyped_n; j++){
			p1 = freq_all[nontyped_index_[j]];
			if ( ! Progress::check_abort() ) {
				p.increment();
				if(p1 < maf){
					// impz[nontyped_index_[j]] = 0;
					// r2pred[nontyped_index_[j]] = 0;
				}else{
					for(i = 0; i < m; i++){
						x11 = 0; x12 = 0; x21 = 0; x22 = 0;
						for(k = 0; k < ind; k++){
							gi = bigm[index_[i]][k];
							gj = bigm[nontyped_index_[j]][k];
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
						p1 = x11 + x12;
						p2 = x21 + x22;
						q1 = x11 + x21;
						q2 = x12 + x22;
						r = (x11 * x22 - x12 * x21) / (sqrt(p1 * p2 * q1 * q2));
						Rit[i] = r;
					}
					beta = Rit.t() * iRtt;
					impz.row(nontyped_index_[j]) = beta * eff_typed;
				}
			}
		}
	}else{

		#pragma omp parallel for schedule(dynamic) firstprivate(Rit) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r, beta)
		for(j = 0; j < nontyped_n; j++){

			p1 = freq_all[nontyped_index_[j]];
			m1 = mean_all[nontyped_index_[j]];
			s1 = sum_all[nontyped_index_[j]];

			if ( ! Progress::check_abort() ) {

				p.increment();

				if(p1 < maf){

					// not impute
					// impz[nontyped_index_[j]] = 0;
					// r2pred[nontyped_index_[j]] = 0;

				}else{

					// impute
					p1 = sd_all[nontyped_index_[j]];
					for(i = 0; i < m; i++){
						p12 = 0;
						p2 = sd_all[index_[i]];
						m2 = mean_all[index_[i]];
						s2 = sum_all[index_[i]];
						for(k = 0; k < ind; k++){
							p12 += (bigm[index_[i]][k]) * (bigm[nontyped_index_[j]][k]);
						}
						p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
						r = p12 / (p1 * p2);
						// if(sparse){
						// 	if(r * r * ind < chisq){
						// 		r = 0.0;
						// 	}
						// }
						Rit[i] = r;
					}

					// get beta
					beta = Rit.t() * iRtt;
					// impute missing and get var
					// impz[nontyped_index_[j], 1] = arma::as_scalar(beta * eff_typed);
					impz.row(nontyped_index_[j]) = beta * eff_typed;
					// r2pred.row(nontyped_index_[j]) = beta * Rit;
				}
			}
		}
	}
	// return List::create(Named("impz") = impz, Named("r2pred") = r2pred);
	return wrap(impz);
}

// [[Rcpp::export]]
SEXP SImputeZ_ld_bin_c(SEXP pBigMat, SEXP pBigMat_G, const IntegerVector typed_index, const arma::mat typed_value, const double lambda = 0.0, const double maf = 0.01, const bool haps = false, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);
	XPtr<BigMatrix> xpMat_G(pBigMat_G);

	switch(xpMat->matrix_type()) {
	case 1:
		return SImputeZ_ld_bin_c<char>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	case 2:
		return SImputeZ_ld_bin_c<short>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	case 4:
		return SImputeZ_ld_bin_c<int>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	case 8:
		return SImputeZ_ld_bin_c<double>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, haps, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

// template <typename T>
// arma::mat SImputeB_Rii_c(XPtr<BigMatrix> pMat, const NumericVector sd_all, const NumericVector mean_all, const NumericVector sum_all, const bool verbose = true, const int threads = 0){
	
// 	omp_setup(threads);

// 	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

// 	int n = pMat->ncol();
// 	int ind = pMat->nrow();
// 	int i, j, k;

// 	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0;
	
// 	arma::mat W(n, n);
// 	MinimalProgressBar pb;
// 	Progress p(n, verbose, pb);

// 	#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r)
// 	for (j = 0; j < n; j++){
// 		if ( ! Progress::check_abort() ) {
// 			p.increment();
// 			p1 = sd_all[j];
// 			m1 = mean_all[j];
// 			s1 = sum_all[j];
// 			W(j, j) = 1;
// 			for(i = (j + 1); i < n; i++){
// 				p12 = 0;
// 				p2 = sd_all[i];
// 				m2 = mean_all[i];
// 				s2 = sum_all[i];
// 				for(k = 0; k < ind; k++){
// 					p12 += (bigm[i][k]) * (bigm[j][k]);
// 				}
// 				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
// 				W(i, j) = p12 / (p2 * p2);
// 				W(j, i) = p12 / (p1 * p1);
// 			}
// 		}
// 	}
// 	return W * W.t();
// }

// // [[Rcpp::export]]
// arma::mat SImputeB_Rii_c(SEXP pBigMat, const NumericVector sd_all, const NumericVector mean_all, const NumericVector sum_all, const bool verbose = true, const int threads = 0){

// 	XPtr<BigMatrix> xpMat(pBigMat);

// 	switch(xpMat->matrix_type()) {
// 	case 1:
// 		return SImputeB_Rii_c<char>(xpMat, sd_all, mean_all, sum_all, verbose, threads);
// 	case 2:
// 		return SImputeB_Rii_c<short>(xpMat, sd_all, mean_all, sum_all, verbose, threads);
// 	case 4:
// 		return SImputeB_Rii_c<int>(xpMat, sd_all, mean_all, sum_all, verbose, threads);
// 	case 8:
// 		return SImputeB_Rii_c<double>(xpMat, sd_all, mean_all, sum_all, verbose, threads);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

// template <typename T>
// SEXP SImputeB_bin_c(XPtr<BigMatrix> pMat, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	omp_setup(threads);

// 	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

// 	int n = pMat->ncol();
// 	int ind = pMat->nrow();
// 	int m, i, j, k;

// 	if(lambda.length() != typed_value.n_cols){
// 		throw Rcpp::exception("length of lambda not equals to Number of traits!");
// 	}

// 	IntegerVector index = seq(0, n - 1);
// 	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
// 	// arma::mat r2pred(n, typed_value.n_cols);
// 	IntegerVector typed_index_ = typed_index - 1;
// 	IntegerVector nontyped_index_ = rev(Rcpp::setdiff(index, typed_index_));

// 	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index_.size() << endl;
// 	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

// 	if(verbose)	Rcerr << "MAF calculating" << endl;
// 	NumericVector freq_all = freq(pMat, index, threads);

// 	// // for hap
// 	// NumericVector freq_all = freq_hap(pMat, index, threads);

// 	List Stat = BigStat(pMat, index, threads);
// 	NumericVector mean_all = Stat[0];
// 	NumericVector sd_all = Stat[1];
// 	NumericVector sum_all = Stat[2];

// 	IntegerVector eff_index = which_c(freq_all[typed_index_], maf, 6);

// 	for(i = 0; i < typed_index_.size(); i++){
// 		for(j = 0; j < typed_value.n_cols; j++){
// 			impz(typed_index_[i], j) = typed_value(i, j);
// 			// r2pred[typed_index_[i], j] = 1;
// 		}
// 	}

// 	IntegerVector index_ = typed_index_[eff_index];
// 	arma::mat eff_typed = typed_value.rows(as<uvec>(eff_index));
// 	m = index_.size();

// 	// // for hap
// 	// NumericVector typed_freq = freq_all[index_];

// 	if(verbose)	Rcerr << "Constructing Var-Cov matrix for all " << n << " SNPs" << endl;
// 	arma::mat Rii = SImputeB_Rii_c(pMat, sd_all, mean_all, sum_all, verbose, threads);
// 	arma::mat Rit = Rii.submat(as<uvec>(nontyped_index_), as<uvec>(index_));

// 	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;
// 	if(verbose)	Rcerr << "Eigen decomposition of sigma_tt matrix" << endl;
// 	arma::vec eigval;
// 	arma::mat eigvec;
// 	eig_sym(eigval, eigvec, Rii.submat(as<uvec>(index_), as<uvec>(index_)));
	
// 	if(verbose)	Rcerr << "Imputing started..." << endl;
// 	for(i = 0; i < eff_typed.n_cols; i++){
// 		arma::rowvec eigveclambda = (1 / (eigval + lambda[i])).t();
// 		arma::mat eiglam = eigvec.each_row() % eigveclambda;
// 		arma::mat iRtt = eiglam * eigvec.t();
// 		arma::vec ipt = Rit * iRtt * eff_typed.col(i);

// 		for(j = 0; j < nontyped_index_.size(); j++){
// 			double p1 = freq_all[nontyped_index_[j]];
// 			if(p1 < maf || p1 > (1 - maf)){

// 				// not impute
// 				// impz[nontyped_index_[j]] = 0;
// 				// r2pred[nontyped_index_[j]] = 0;
// 			}else{

// 				// impute
// 				impz(nontyped_index_[j], i) = ipt[j];
// 				// r2pred.row(nontyped_index_[j]) = beta * Rit;
// 			}
// 		}
// 	}

// 	// return List::create(Named("impz") = impz, Named("r2pred") = r2pred);
// 	return wrap(impz);
// }

// // [[Rcpp::export]]
// SEXP SImputeB_bin_c(SEXP pBigMat, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	XPtr<BigMatrix> xpMat(pBigMat);

// 	switch(xpMat->matrix_type()) {
// 	case 1:
// 		return SImputeB_bin_c<char>(xpMat, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 2:
// 		return SImputeB_bin_c<short>(xpMat, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 4:
// 		return SImputeB_bin_c<int>(xpMat, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 8:
// 		return SImputeB_bin_c<double>(xpMat, typed_index, typed_value, lambda, maf, verbose, threads);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

// template <typename T>
// arma::mat SImputeB_ld_Rii_c(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pMat_G, const LogicalVector refindx, const IntegerVector typed_index_, const NumericVector sd_all, const NumericVector mean_all, const NumericVector sum_all, const NumericVector sd_gwas, const NumericVector mean_gwas, const NumericVector sum_gwas, const bool verbose = true, const int threads = 0){
	
// 	omp_setup(threads);

// 	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
// 	MatrixAccessor<T> bigm_G = MatrixAccessor<T>(*pMat_G);

// 	int n = pMat->ncol();
// 	int ind = pMat->nrow();
// 	int ngwas = pMat_G->ncol();
// 	int indgwas = pMat_G->nrow();
// 	int i, j, k;

// 	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0;
	
// 	arma::mat W(n, n);

// 	for(int xx = 0; xx < 2; xx++){

// 		if(xx == 0){

// 			MinimalProgressBar pb;
// 			Progress p(n, verbose, pb);

// 			#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r)
// 			for (j = 0; j < n; j++){
// 				if ( ! Progress::check_abort() ) {
// 					p.increment();
// 					p1 = sd_all[j];
// 					m1 = mean_all[j];
// 					s1 = sum_all[j];
// 					W(j, j) = 1;
// 					for(i = (j + 1); i < n; i++){
// 						if(refindx[j] && refindx[i]){
// 								// nothing
// 						}else{
// 							p12 = 0;
// 							p2 = sd_all[i];
// 							m2 = mean_all[i];
// 							s2 = sum_all[i];
// 							for(k = 0; k < ind; k++){
// 								p12 += (bigm[i][k]) * (bigm[j][k]);
// 							}
// 							p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
// 							W(i, j) = p12 / (p2 * p2);
// 							W(j, i) = p12 / (p1 * p1);
// 						}
// 					}
// 				}
// 			}
// 		}else{

// 			if(verbose)	Rcerr << "Update W for typed SNPs" << endl;
// 			MinimalProgressBar pb;
// 			Progress p(ngwas, verbose, pb);

// 			#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r)
// 			for (j = 0; j < ngwas; j++){
// 				if ( ! Progress::check_abort() ) {
// 					p.increment();
// 					p1 = sd_gwas[j];
// 					m1 = mean_gwas[j];
// 					s1 = sum_gwas[j];
// 					for(i = (j + 1); i < ngwas; i++){
// 						p12 = 0;
// 						p2 = sd_gwas[i];
// 						m2 = mean_gwas[i];
// 						s2 = sum_gwas[i];
// 						for(k = 0; k < indgwas; k++){
// 							p12 += (bigm_G[i][k]) * (bigm_G[j][k]);
// 						}
// 						p12 -= s1 * m2 + s2 * m1 - indgwas * m1 * m2;
// 						W(typed_index_[i], typed_index_[j]) = p12 / (p2 * p2);
// 						W(typed_index_[j], typed_index_[i]) = p12 / (p1 * p1);
// 					}
// 				}
// 			}
// 		}
// 	}
// 	return W * W.t();
// }

// // [[Rcpp::export]]
// arma::mat SImputeB_ld_Rii_c(SEXP pBigMat, SEXP pBigMat_G, const LogicalVector refindx, const IntegerVector typed_index_, const NumericVector sd_all, const NumericVector mean_all, const NumericVector sum_all, const NumericVector sd_gwas, const NumericVector mean_gwas, const NumericVector sum_gwas, const bool verbose = true, const int threads = 0){

// 	XPtr<BigMatrix> xpMat(pBigMat);
// 	XPtr<BigMatrix> xpMat_G(pBigMat_G);

// 	switch(xpMat->matrix_type()) {
// 	case 1:
// 		return SImputeB_ld_Rii_c<char>(xpMat, xpMat_G, refindx, typed_index_, sd_all, mean_all, sum_all, sd_gwas, mean_gwas, sum_gwas, verbose, threads);
// 	case 2:
// 		return SImputeB_ld_Rii_c<short>(xpMat, xpMat_G, refindx, typed_index_, sd_all, mean_all, sum_all, sd_gwas, mean_gwas, sum_gwas, verbose, threads);
// 	case 4:
// 		return SImputeB_ld_Rii_c<int>(xpMat, xpMat_G, refindx, typed_index_, sd_all, mean_all, sum_all, sd_gwas, mean_gwas, sum_gwas, verbose, threads);
// 	case 8:
// 		return SImputeB_ld_Rii_c<double>(xpMat, xpMat_G, refindx, typed_index_, sd_all, mean_all, sum_all, sd_gwas, mean_gwas, sum_gwas, verbose, threads);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

// template <typename T>
// SEXP SImputeB_ld_bin_c(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pMat_G, const LogicalVector refindx, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	omp_setup(threads);

// 	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
// 	MatrixAccessor<T> bigm_G = MatrixAccessor<T>(*pMat_G);

// 	if(lambda.length() != typed_value.n_cols){
// 		throw Rcpp::exception("length of lambda not equals to Number of traits!");
// 	}

// 	int n = pMat->ncol();
// 	int n_gwas = pMat_G->ncol();
// 	int ind = pMat->nrow();
// 	int ind_gwas = pMat_G->nrow();
// 	int m, i, j, k;
// 	double p1 = 0.0, m1 = 0.0, s1 = 0.1, p2 = 0.0, m2 = 0.0, s2 = 0, p12 = 0.0, r = 0.0;
// 	// arma::mat typed_value = Rcpp::as<arma::mat>(typed_value_R);
// 	IntegerVector index = seq(0, n - 1);
// 	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
// 	// arma::mat r2pred(n, typed_value.n_cols);
// 	MinimalProgressBar pb;
// 	IntegerVector typed_index_ = typed_index - 1;
// 	IntegerVector nontyped_index_ = rev(Rcpp::setdiff(index, typed_index_));

// 	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index_.size() << endl;
// 	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

// 	if(verbose)	Rcerr << "MAF calculating" << endl;
// 	NumericVector freq_all = freq(pMat, index, threads);
// 	// NumericVector freq_gwas_all = freq(pMat_G, seq(0, n_gwas - 1), threads);
// 	List Stat = BigStat(pMat, index, threads);
// 	List Stat_gwas = BigStat(pMat_G, seq(0, n_gwas - 1), threads);
// 	NumericVector mean_all = Stat[0];
// 	NumericVector sd_all = Stat[1];
// 	NumericVector sum_all = Stat[2];
// 	NumericVector mean_gwas = Stat_gwas[0];
// 	NumericVector sd_gwas = Stat_gwas[1];
// 	NumericVector sum_gwas = Stat_gwas[2];

// 	IntegerVector eff_index = which_c(freq_all[typed_index_], maf, 6);

// 	// typed SNPs needn't to impute
// 	for(i = 0; i < typed_index_.size(); i++){
// 		for(j = 0; j < typed_value.n_cols; j++){
// 			impz(typed_index_[i], j) = typed_value(i, j);
// 			// r2pred[typed_index_[i], j] = 1;
// 		}
// 	}

// 	IntegerVector index_ = typed_index_[eff_index];
// 	arma::mat eff_typed(eff_index.size(), typed_value.n_cols);
// 	for(i = 0; i < eff_index.size(); i++){
// 		eff_typed.row(i) = typed_value.row(eff_index[i]);
// 	}
// 	m = index_.size();

// 	if(verbose)	Rcerr << "Constructing Var-Cov matrix for all " << n << " SNPs" << endl;
// 	arma::mat Rii = SImputeB_ld_Rii_c(pMat, pMat_G, refindx, typed_index_, sd_all, mean_all, sum_all, sd_gwas, mean_gwas, sum_gwas, verbose, threads);
// 	arma::mat Rit = Rii.submat(as<uvec>(nontyped_index_), as<uvec>(index_));

// 	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;
// 	if(verbose)	Rcerr << "Eigen decomposition of sigma_tt matrix" << endl;
// 	arma::vec eigval;
// 	arma::mat eigvec;
// 	eig_sym(eigval, eigvec, Rii.submat(as<uvec>(index_), as<uvec>(index_)));
	
// 	if(verbose)	Rcerr << "Imputing started..." << endl;
// 	for(i = 0; i < eff_typed.n_cols; i++){
// 		arma::rowvec eigveclambda = (1 / (eigval + lambda[i])).t();
// 		arma::mat eiglam = eigvec.each_row() % eigveclambda;
// 		arma::mat iRtt = eiglam * eigvec.t();
// 		arma::vec ipt = Rit * iRtt * eff_typed.col(i);

// 		for(j = 0; j < nontyped_index_.size(); j++){
// 			double p1 = freq_all[nontyped_index_[j]];
// 			if(p1 < maf || p1 > (1 - maf)){

// 				// not impute
// 				// impz[nontyped_index_[j]] = 0;
// 				// r2pred[nontyped_index_[j]] = 0;
// 			}else{

// 				// impute
// 				impz(nontyped_index_[j], i) = ipt[j];
// 				// r2pred.row(nontyped_index_[j]) = beta * Rit;
// 			}
// 		}
// 	}

// 	return wrap(impz);
// }

// // [[Rcpp::export]]
// SEXP SImputeB_ld_bin_c(SEXP pBigMat, SEXP pBigMat_G, const LogicalVector refindx, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	XPtr<BigMatrix> xpMat(pBigMat);
// 	XPtr<BigMatrix> xpMat_G(pBigMat_G);

// 	switch(xpMat->matrix_type()) {
// 	case 1:
// 		return SImputeB_ld_bin_c<char>(xpMat, xpMat_G, refindx, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 2:
// 		return SImputeB_ld_bin_c<short>(xpMat, xpMat_G, refindx, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 4:
// 		return SImputeB_ld_bin_c<int>(xpMat, xpMat_G, refindx, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 8:
// 		return SImputeB_ld_bin_c<double>(xpMat, xpMat_G, refindx, typed_index, typed_value, lambda, maf, verbose, threads);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

// template <typename T>
// SEXP SImputeB_new_bin_c(XPtr<BigMatrix> pMat, const int n_gwas, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	omp_setup(threads);

// 	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

// 	if(lambda.length() != typed_value.n_cols){
// 		throw Rcpp::exception("length of lambda not equals to Number of traits!");
// 	}

// 	int n = pMat->ncol();
// 	int ind = pMat->nrow();
// 	int m, i, j, k;
// 	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0;
// 	IntegerVector index = seq(0, n - 1);
// 	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
// 	MinimalProgressBar pb;
// 	IntegerVector typed_index_ = typed_index - 1;
// 	IntegerVector nontyped_index_ = rev(Rcpp::setdiff(index, typed_index_));

// 	// typed SNPs needn't to impute
// 	for(i = 0; i < typed_index_.size(); i++){
// 		for(j = 0; j < typed_value.n_cols; j++){
// 			impz(typed_index_[i], j) = typed_value(i, j);
// 		}
// 	}

// 	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index_.size() << endl;
// 	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

// 	if(verbose)	Rcerr << "MAF calculating" << endl;
// 	NumericVector freq_all = freq(pMat, index, threads);
// 	List Stat = BigStat(pMat, index, threads);
// 	NumericVector mean_all = Stat[0];
// 	NumericVector sd_all = Stat[1];
// 	NumericVector sum_all = Stat[2];

// 	IntegerVector eff_index = which_c(freq_all[typed_index_], maf, 6);
// 	IntegerVector index_ = typed_index_[eff_index];
// 	NumericVector typed_sd = sd_all[index_];
// 	NumericVector typed_mean = mean_all[index_];
// 	NumericVector typed_sum = sum_all[index_];
// 	arma::mat eff_typed = typed_value.rows(as<uvec>(eff_index));
// 	m = index_.size();

// 	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;

// 	arma::mat Rtt(m, m);

// 	#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r)
// 	for (j = 0; j < m; j++){
// 		p1 = typed_sd[j];
// 		m1 = typed_mean[j];
// 		s1 = typed_sum[j];
// 		for(i = j + 1; i < m; i++){
// 			p12 = 0;
// 			p2 = typed_sd[i];
// 			m2 = typed_mean[i];
// 			s2 = typed_sum[i];
// 			for(k = 0; k < ind; k++){
// 				p12 += (bigm[index_[i]][k]) * (bigm[index_[j]][k]);
// 			}
// 			p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
// 			// r = p12 / (p1 * p2);
// 			Rtt(i, j) = Rtt(j, i) = p12 * n_gwas / ind;
// 		}
// 		Rtt(j, j) = p1 * p1 * n_gwas / ind;
// 	}

// 	if(verbose)	Rcerr << "Eigen decomposition of sigma_tt matrix" << endl;
// 	arma::vec eigval;
// 	arma::mat eigvec;
// 	eig_sym(eigval, eigvec, Rtt);
	
// 	if(verbose)	Rcerr << "Estimate joint effect" << endl;
// 	for(int i = 0; i < eff_typed.n_rows; i++){
// 		eff_typed.row(i) *= Rtt(i, i);
// 	}

// 	arma::mat joint_eff(eff_typed.n_rows, eff_typed.n_cols);
// 	for(i = 0; i < eff_typed.n_cols; i++){
// 		arma::rowvec eigveclambda = (1 / (eigval + lambda[i])).t();
// 		arma::mat eiglam = eigvec.each_row() % eigveclambda;
// 		arma::mat iRtt = eiglam * eigvec.t();
// 		joint_eff.col(i) = iRtt * eff_typed.col(i);
// 	}

// 	// string outfile1 = "snp.blup";
// 	// ofstream out1;
// 	// out1.open(outfile1, ios::out|ios::app);
// 	// for(int i = 0; i < eff_typed.n_rows; i++){
// 	// 	out1 << boost::format("%= 12.4f %= 12.4f %= 12.4f\n")
// 	// 	% joint_eff(i, 0)
// 	//     % joint_eff(i, 1)
// 	//     % joint_eff(i, 2);
// 	// }
// 	// out1.close();

// 	// arma::vec Rit(m);
// 	if(verbose)	Rcerr << "Imputing started..." << endl;

// 	Progress p(nontyped_index_.size(), verbose, pb);

// 	#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, i, k, p12, p2, m2, s2, r)
// 	for(j = 0; j < nontyped_index_.size(); j++){
// 		if(freq_all[nontyped_index_[j]] < maf || freq_all[nontyped_index_[j]] > (1 - maf))	continue;
// 		arma::vec Rit(m);
// 		p1 = sd_all[nontyped_index_[j]];
// 		m1 = mean_all[nontyped_index_[j]];
// 		s1 = sum_all[nontyped_index_[j]];
// 		if ( ! Progress::check_abort() ) {
// 			p.increment();
// 			for(i = 0; i < m; i++){
// 				p12 = 0;
// 				p2 = sd_all[index_[i]];
// 				m2 = mean_all[index_[i]];
// 				s2 = sum_all[index_[i]];
// 				for(k = 0; k < ind; k++){
// 					p12 += (bigm[index_[i]][k]) * (bigm[nontyped_index_[j]][k]);
// 				}
// 				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
// 				r = p12 / (p1 * p2);
// 				Rit[i] = (1 / p1) * r * p2;
// 			}
// 			impz.row(nontyped_index_[j]) = Rit.t() * joint_eff;
// 		}
// 	}
// 	return wrap(impz);
// }

// // [[Rcpp::export]]
// SEXP SImputeB_new_bin_c(SEXP pBigMat, const int n_gwas, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	XPtr<BigMatrix> xpMat(pBigMat);

// 	switch(xpMat->matrix_type()) {
// 	case 1:
// 		return SImputeB_new_bin_c<char>(xpMat, n_gwas, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 2:
// 		return SImputeB_new_bin_c<short>(xpMat, n_gwas, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 4:
// 		return SImputeB_new_bin_c<int>(xpMat, n_gwas, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 8:
// 		return SImputeB_new_bin_c<double>(xpMat, n_gwas, typed_index, typed_value, lambda, maf, verbose, threads);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

// template <typename T>
// SEXP SImputeB_ld_new_bin_c(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pMat_G, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	omp_setup(threads);

// 	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
// 	MatrixAccessor<T> bigm_G = MatrixAccessor<T>(*pMat_G);

// 	if(lambda.length() != typed_value.n_cols){
// 		throw Rcpp::exception("length of lambda not equals to Number of traits!");
// 	}

// 	int n = pMat->ncol();
// 	int ind = pMat->nrow();
// 	int ngwas = pMat_G->ncol();
// 	int indgwas = pMat_G->nrow();
// 	int m, i, j, k;
// 	double p1 = 0.0, m1 = 0.0, s1 = 0.0, p2 = 0.0, m2 = 0.0, s2 = 0.0, p12 = 0.0, r = 0.0;
// 	IntegerVector index = seq(0, n - 1);
// 	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
// 	MinimalProgressBar pb;
// 	IntegerVector typed_index_ = typed_index - 1;
// 	IntegerVector nontyped_index_ = rev(Rcpp::setdiff(index, typed_index_));

// 	// typed SNPs needn't to impute
// 	for(i = 0; i < typed_index_.size(); i++){
// 		for(j = 0; j < typed_value.n_cols; j++){
// 			impz(typed_index_[i], j) = typed_value(i, j);
// 		}
// 	}
	
// 	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index_.size() << endl;
// 	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

// 	if(verbose)	Rcerr << "MAF calculating" << endl;
// 	NumericVector freq_all = freq(pMat, index, threads);
// 	List Stat = BigStat(pMat, index, threads);
// 	List Stat_gwas = BigStat(pMat_G, seq(0, ngwas - 1), threads);
// 	NumericVector mean_all = Stat[0];
// 	NumericVector sd_all = Stat[1];
// 	NumericVector sum_all = Stat[2];
// 	NumericVector mean_gwas = Stat_gwas[0];
// 	NumericVector sd_gwas = Stat_gwas[1];
// 	NumericVector sum_gwas = Stat_gwas[2];

// 	IntegerVector eff_index = which_c(freq_all[typed_index_], maf, 6);
// 	IntegerVector index_ = typed_index_[eff_index];
// 	NumericVector typed_sd = sd_gwas[eff_index];
// 	NumericVector typed_mean = mean_gwas[eff_index];
// 	NumericVector typed_sum = sum_gwas[eff_index];
// 	arma::mat eff_typed = typed_value.rows(as<uvec>(eff_index));
// 	m = index_.size();

// 	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;

// 	arma::mat Rtt(m, m);

// 	#pragma omp parallel for schedule(dynamic) private(j, s1, s2, p1, m1, i, k, p12, p2, m2, r)
// 	for (j = 0; j < m; j++){
// 		p1 = typed_sd[j];
// 		m1 = typed_mean[j];
// 		s1 = typed_sum[j];
// 		for(i = j + 1; i < m; i++){
// 			p12 = 0;
// 			p2 = typed_sd[i];
// 			m2 = typed_mean[i];
// 			s2 = typed_sum[i];
// 			for(k = 0; k < indgwas; k++){
// 				p12 += (bigm_G[eff_index[i]][k]) * (bigm_G[eff_index[j]][k]);
// 			}
// 			p12 -= s1 * m2 + s2 * m1 - indgwas * m1 * m2;
// 			// r = p12 / (p1 * p2);
// 			Rtt(i, j) = Rtt(j, i) = p12;
// 		}
// 		Rtt(j, j) = p1 * p1;
// 	}

// 	if(verbose)	Rcerr << "Eigen decomposition of sigma_tt matrix" << endl;
// 	arma::vec eigval;
// 	arma::mat eigvec;
// 	eig_sym(eigval, eigvec, Rtt);
	
// 	if(verbose)	Rcerr << "Estimate joint effect" << endl;
// 	for(int i = 0; i < eff_typed.n_rows; i++){
// 		eff_typed.row(i) *= Rtt(i, i);
// 	}

// 	arma::mat joint_eff(eff_typed.n_rows, eff_typed.n_cols);
// 	for(i = 0; i < eff_typed.n_cols; i++){
// 		arma::rowvec eigveclambda = (1 / (eigval + lambda[i])).t();
// 		arma::mat eiglam = eigvec.each_row() % eigveclambda;
// 		arma::mat iRtt = eiglam * eigvec.t();
// 		joint_eff.col(i) = iRtt * eff_typed.col(i);
// 	}

// 	// string outfile1 = "snp.blup";
// 	// ofstream out1;
// 	// out1.open(outfile1, ios::out|ios::app);
// 	// for(int i = 0; i < eff_typed.n_rows; i++){
// 	// 	out1 << boost::format("%= 12.4f %= 12.4f %= 12.4f\n")
// 	// 	% joint_eff(i, 0)
// 	//     % joint_eff(i, 1)
// 	//     % joint_eff(i, 2);
// 	// }
// 	// out1.close();

// 	// arma::vec Rit(m);
// 	if(verbose)	Rcerr << "Imputing started..." << endl;

// 	Progress p(nontyped_index_.size(), verbose, pb);

// 	#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, i, k, p12, p2, m2, s2, r)
// 	for(j = 0; j < nontyped_index_.size(); j++){
// 		if(freq_all[nontyped_index_[j]] < maf || freq_all[nontyped_index_[j]] > (1 - maf))	continue;
// 		arma::vec Rit(m);
// 		p1 = sd_all[nontyped_index_[j]];
// 		m1 = mean_all[nontyped_index_[j]];
// 		s1 = sum_all[nontyped_index_[j]];
// 		if ( ! Progress::check_abort() ) {
// 			p.increment();
// 			for(i = 0; i < m; i++){
// 				p12 = 0;
// 				p2 = sd_all[index_[i]];
// 				m2 = mean_all[index_[i]];
// 				s2 = sum_all[index_[i]];
// 				for(k = 0; k < ind; k++){
// 					p12 += (bigm[index_[i]][k]) * (bigm[nontyped_index_[j]][k]);
// 				}
// 				p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
// 				r = p12 / (p1 * p2);
// 				Rit[i] = (1 / p1) * r * p2;
// 			}
// 			impz.row(nontyped_index_[j]) = Rit.t() * joint_eff;
// 		}
// 	}
// 	return wrap(impz);
// }

// // [[Rcpp::export]]
// SEXP SImputeB_ld_new_bin_c(SEXP pBigMat, SEXP pBigMat_G, const IntegerVector typed_index, const arma::mat typed_value, const NumericVector lambda, const double maf = 0.01, const bool verbose = true, const int threads = 0){

// 	XPtr<BigMatrix> xpMat(pBigMat);
// 	XPtr<BigMatrix> xpMat_G(pBigMat_G);

// 	switch(xpMat->matrix_type()) {
// 	case 1:
// 		return SImputeB_ld_new_bin_c<char>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 2:
// 		return SImputeB_ld_new_bin_c<short>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 4:
// 		return SImputeB_ld_new_bin_c<int>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, verbose, threads);
// 	case 8:
// 		return SImputeB_ld_new_bin_c<double>(xpMat, xpMat_G, typed_index, typed_value, lambda, maf, verbose, threads);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

template <typename T>
SEXP SImpute_joint_c(XPtr<BigMatrix> pMat, const IntegerVector typed_index, const arma::mat typed_value, const double maf = 0.01, const bool verbose = true, int threads = 0){

	omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int n = pMat->ncol();
	int ind = pMat->nrow();
	int m = typed_index.size(), i, j, k;
	double p1 = 0.0, m1, m2, p2 = 0.0, p12 = 0.0, r = 0.0, s1, s2;
	MinimalProgressBar pb;
	// MinimalProgressBar_plus pb;
	IntegerVector typed_index_ = typed_index - 1;
	IntegerVector index = seq(0, n - 1);
	// arma::mat typed_value = Rcpp::as<arma::mat>(typed_value_R);
	arma::mat impz(n, typed_value.n_cols); impz.fill(0);

	if(verbose)	Rcerr << "Number of typed SNPs is " << typed_index_.size() << endl;
	if(verbose)	Rcerr << "Number of Total SNPs is " << n << endl;

	if(verbose)	Rcerr << "MAF calculating" << endl;
	NumericVector freq_all = freq(pMat, index, threads = threads);
	List Stat = BigStat(pMat, index, threads = threads);
	NumericVector mean_all = Stat[0];
	NumericVector sd_all = Stat[1];
	NumericVector sum_all = Stat[2];

	if(verbose)	Rcerr << "Imputing started..." << endl;

	Progress p(n, verbose, pb);

	#pragma omp parallel for schedule(dynamic) private(j, p1, m1, s1, i, k, p12, p2, m2, s2, r)
	for(j = 0; j < n; j++){
		if(freq_all[j] < maf)	continue;
		arma::vec Rit(m);
		p1 = sd_all[j];
		m1 = mean_all[j];
		s1 = sum_all[j];
		if ( ! Progress::check_abort() ) {
			p.increment();
			for(i = 0; i < m; i++){
				if(typed_index_[i] == j){
					Rit[i] = 1;
				}else{
					p12 = 0;
					p2 = sd_all[typed_index_[i]];
					m2 = mean_all[typed_index_[i]];
					s2 = sum_all[typed_index_[i]];
					for(k = 0; k < ind; k++){
						p12 += (bigm[typed_index_[i]][k]) * (bigm[j][k]);
					}
					p12 -= s1 * m2 + s2 * m1 - ind * m1 * m2;
					r = p12 / (p1 * p2);
					Rit[i] = (1 / p1) * r * p2;
				}
			}
			impz.row(j) = Rit.t() * typed_value;
		}
	}
	return wrap(impz);
}

// [[Rcpp::export]]
SEXP SImpute_joint_c(SEXP pBigMat, const IntegerVector typed_index, const arma::mat typed_value, const double maf = 0.01, const bool verbose = true, int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return SImpute_joint_c<char>(xpMat, typed_index, typed_value, maf, verbose, threads);
	case 2:
		return SImpute_joint_c<short>(xpMat, typed_index, typed_value, maf, verbose, threads);
	case 4:
		return SImpute_joint_c<int>(xpMat, typed_index, typed_value, maf, verbose, threads);
	case 8:
		return SImpute_joint_c<double>(xpMat, typed_index, typed_value, maf, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
