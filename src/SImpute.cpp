#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#include "stats.h"
#include "probar.h"

template <typename T>
SEXP SImputeZ_bin_c(XPtr<BigMatrix> pMat, const IntegerVector & bin_index, const IntegerVector & typed_index, const arma::mat & typed_value, const double lambda = 0.0, const double maf = 0.01, const bool verbose = true, const int threads = 0){

	omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int n = bin_index.size();
	int ind = pMat->nrow();
	int m, i, j, k, g1, g2;

	double p1 = 0.0, p2 = 0.0, p12 = 0.0, r = 0.0;
	// arma::mat typed_value = Rcpp::as<arma::mat>(typed_value_R);
	IntegerVector index = seq(0, n - 1);
	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
	// arma::mat r2pred(n, typed_value.n_cols);
	MinimalProgressBar pb;
	// MinimalProgressBar_plus pb;
	// IntegerVector typed_index_ = typed_index - 1;
	IntegerVector nontyped_index = rev(Rcpp::setdiff(index, typed_index));
	IntegerVector bin_nontyped_index = bin_index[nontyped_index];

	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index.size() << endl;
	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

	if(verbose)	Rcerr << "MAF calculating" << endl;
	NumericVector freq_all = freq_snp(pMat, bin_index, threads);

	IntegerVector eff_index = which_c(freq_all[typed_index], maf, 6);

	int typed_n = typed_index.size();
	int typed_ncol = typed_value.n_cols;
	for(i = 0; i < typed_n; i++){
		for(j = 0; j < typed_ncol; j++){
			impz(typed_index[i], j) = typed_value(i, j);
			// r2pred[typed_index_[i], j] = 1;
		}
	}

	IntegerVector typed_eff_index = typed_index[eff_index];
	IntegerVector bin_typed_index = bin_index[typed_eff_index];
	arma::mat eff_typed = typed_value.rows(as<uvec>(eff_index));
	m = bin_typed_index.size();

	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;

	arma::mat Rtt(m, m);

	#pragma omp parallel for schedule(dynamic) private(j, i, k, p1, p2, g1, g2, p12, r)
	for (j = 0; j < m; j++){
		p1 = freq_all[typed_eff_index[j]];
		for(i = j + 1; i < m; i++){
			p2 = freq_all[typed_eff_index[i]];
			p12 = 0.0;
			double fm = sqrt(p1 * (1.0 - p1) * p2 * (1.0 - p2));
			if(fm){
				for(k = 0; k < ind; k++){
					g1 = bigm[bin_typed_index[i]][k];
					g2 = bigm[bin_typed_index[j]][k];
					p12 += (g1 == 11 && g2 == 11) ? (2) : ((g1 == 22 || g2 == 22 || ((g1 + g2) == 33)) ? (0) : (1));
				}
				p12 /= (ind * 2);
				r = (p12 - p1 * p2) / fm;
			}else{
				r = 0.0;
			}
			Rtt(i, j) = Rtt(j, i) = r;
		}
		Rtt(j, j) = 1 + lambda;
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
	int nontyped_n = nontyped_index.size();
	Progress p(nontyped_n, verbose, pb);
	arma::vec Rit(m); Rit.fill(0);

	#pragma omp parallel for schedule(dynamic) firstprivate(Rit) private(j, i, k, p1, p2, p12, g1, g2, r, beta)
	for(j = 0; j < nontyped_n; j++){
		p1 = freq_all[nontyped_index[j]];
		if ( ! Progress::check_abort() ) {
			p.increment();
			if(p1 < maf || p1 > (1 - maf)){
				// impz[nontyped_index_[j]] = 0;
				// r2pred[nontyped_index_[j]] = 0;
			}else{
				for(i = 0; i < m; i++){
					p2 = freq_all[typed_eff_index[i]];
					p12 = 0.0;
					double fm = sqrt(p1 * (1.0 - p1) * p2 * (1.0 - p2));
					if(fm){
						for(k = 0; k < ind; k++){
							g1 = bigm[bin_typed_index[i]][k];
							g2 = bigm[bin_nontyped_index[j]][k];
							p12 += (g1 == 11 && g2 == 11) ? (2) : ((g1 == 22 || g2 == 22 || ((g1 + g2) == 33)) ? (0) : (1));
						}
						p12 /= (ind * 2);
						r = (p12 - p1 * p2) / fm;
					}else{
						r = 0.0;
					}
					Rit[i] = r;
				}
				beta = Rit.t() * iRtt;
				impz.row(nontyped_index[j]) = beta * eff_typed;
			}
		}
	}

	// return List::create(Named("impz") = impz, Named("r2pred") = r2pred);
	return wrap(impz);
}

// [[Rcpp::export]]
SEXP SImputeZ_bin_c(SEXP pBigMat, const IntegerVector & bin_index, const IntegerVector & typed_index, const arma::mat & typed_value, const double lambda = 0.0, const double maf = 0.01, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return SImputeZ_bin_c<char>(xpMat, bin_index, typed_index, typed_value, lambda, maf, verbose, threads);
	case 2:
		return SImputeZ_bin_c<short>(xpMat, bin_index, typed_index, typed_value, lambda, maf, verbose, threads);
	case 4:
		return SImputeZ_bin_c<int>(xpMat, bin_index, typed_index, typed_value, lambda, maf, verbose, threads);
	case 8:
		return SImputeZ_bin_c<double>(xpMat, bin_index, typed_index, typed_value, lambda, maf, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP SImputeZ_ld_bin_c(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pMat_G, const IntegerVector & bin_index , const IntegerVector & typed_index, const IntegerVector & typed_bin_index, const arma::mat & typed_value, const double lambda = 0.0, const double maf = 0.01, const bool verbose = true, const int threads = 0){

	omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);
	MatrixAccessor<T> bigm_G = MatrixAccessor<T>(*pMat_G);

	int n = bin_index.size();
	int ind = pMat->nrow();
	int ind_gwas = pMat_G->nrow();
	int m, i, j, k, g1, g2;
	double p1 = 0.0, p2 = 0.0, p12 = 0.0, r = 0.0;
	// arma::mat typed_value = Rcpp::as<arma::mat>(typed_value_R);
	IntegerVector index = seq(0, n - 1);
	arma::mat impz(n, typed_value.n_cols); impz.fill(0);
	// arma::mat r2pred(n, typed_value.n_cols);
	MinimalProgressBar pb;
	// MinimalProgressBar_plus pb;
	// IntegerVector typed_index_ = typed_index - 1;
	IntegerVector nontyped_index = rev(Rcpp::setdiff(index, typed_index));
	IntegerVector bin_nontyped_index = bin_index[nontyped_index];

	if(verbose)	Rcerr << "Number of typed SNPs of the window is " << typed_index.size() << endl;
	if(verbose)	Rcerr << "Number of Total SNPs of the window is " << n << endl;

	if(verbose)	Rcerr << "MAF calculating" << endl;
	NumericVector freq_all = freq_snp(pMat, bin_index, threads);
	NumericVector freq_all_gwas = freq_snp(pMat_G, typed_bin_index, threads);
	IntegerVector eff_index = which_c(freq_all[typed_index], maf, 6);

	// typed SNPs needn't to impute
	int typed_n = typed_index.size();
	int typed_ncol = typed_value.n_cols;
	for(i = 0; i < typed_n; i++){
		for(j = 0; j < typed_ncol; j++){
			impz(typed_index[i], j) = typed_value(i, j);
			// r2pred[typed_index_[i], j] = 1;
		}
	}

	IntegerVector typed_eff_index = typed_index[eff_index];
	IntegerVector bin_typed_index = bin_index[typed_eff_index];
	IntegerVector gwas_bin_typed_index = typed_bin_index[eff_index];
	NumericVector freq_gwas = freq_all_gwas[eff_index];
	int eff_index_n = eff_index.size();
	arma::mat eff_typed(eff_index_n, typed_value.n_cols);
	for(i = 0; i < eff_index_n; i++){
		eff_typed.row(i) = typed_value.row(eff_index[i]);
	}
	m = bin_typed_index.size();

	if(verbose)	Rcerr << "Using " << m << " SNPs for sigma_tt matrix" << endl;

	arma::mat Rtt(m, m);

	#pragma omp parallel for schedule(dynamic) private(j, i, k, p1, p2, g1, g2, p12, r)
	for (j = 0; j < m; j++){
		p1 = freq_gwas[j];
		for(i = j + 1; i < m; i++){
			p2 = freq_gwas[i];
			p12 = 0.0;
			double fm = sqrt(p1 * (1.0 - p1) * p2 * (1.0 - p2));
			if(fm){
				for(k = 0; k < ind_gwas; k++){
					g1 = bigm_G[gwas_bin_typed_index[i]][k];
					g2 = bigm_G[gwas_bin_typed_index[j]][k];
					p12 += (g1 == 11 && g2 == 11) ? (2) : ((g1 == 22 || g2 == 22 || ((g1 + g2) == 33)) ? (0) : (1));
				}
				p12 /= (ind_gwas * 2);
				r = (p12 - p1 * p2) / fm;
			}else{
				r = 0.0;
			}
			Rtt(i, j) = Rtt(j, i) = r;
		}
		Rtt(j, j) = 1 + lambda;
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

	int nontyped_n = nontyped_index.size();
	Progress p(nontyped_n, verbose, pb);

	#pragma omp parallel for schedule(dynamic) firstprivate(Rit) private(j, i, k, p1, p2, p12, g1, g2, r, beta)
	for(j = 0; j < nontyped_n; j++){
		p1 = freq_all[nontyped_index[j]];
		if ( ! Progress::check_abort() ) {
			p.increment();
			if(p1 < maf || p1 > (1 - maf)){
				// impz[nontyped_index_[j]] = 0;
				// r2pred[nontyped_index_[j]] = 0;
			}else{
				for(i = 0; i < m; i++){
					p2 = freq_all[typed_eff_index[i]];
					p12 = 0.0;
					double fm = sqrt(p1 * (1.0 - p1) * p2 * (1.0 - p2));
					if(fm){
						for(k = 0; k < ind; k++){
							g1 = bigm[bin_typed_index[i]][k];
							g2 = bigm[bin_nontyped_index[j]][k];
							p12 += (g1 == 11 && g2 == 11) ? (2) : ((g1 == 22 || g2 == 22 || ((g1 + g2) == 33)) ? (0) : (1));
						}
						p12 /= (ind * 2);
						r = (p12 - p1 * p2) / fm;
					}else{
						r = 0.0;
					}
					Rit[i] = r;
				}
				beta = Rit.t() * iRtt;
				impz.row(nontyped_index[j]) = beta * eff_typed;
			}
		}
	}

	// return List::create(Named("impz") = impz, Named("r2pred") = r2pred);
	return wrap(impz);
}

// [[Rcpp::export]]
SEXP SImputeZ_ld_bin_c(SEXP pBigMat, SEXP pBigMat_G, const IntegerVector & bin_index , const IntegerVector & typed_index, const IntegerVector & typed_bin_index, const arma::mat & typed_value, const double lambda = 0.0, const double maf = 0.01, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);
	XPtr<BigMatrix> xpMat_G(pBigMat_G);

	switch(xpMat->matrix_type()) {
	case 1:
		return SImputeZ_ld_bin_c<char>(xpMat, xpMat_G, bin_index, typed_index, typed_bin_index, typed_value, lambda, maf, verbose, threads);
	case 2:
		return SImputeZ_ld_bin_c<short>(xpMat, xpMat_G, bin_index, typed_index, typed_bin_index, typed_value, lambda, maf, verbose, threads);
	case 4:
		return SImputeZ_ld_bin_c<int>(xpMat, xpMat_G, bin_index, typed_index, typed_bin_index, typed_value, lambda, maf, verbose, threads);
	case 8:
		return SImputeZ_ld_bin_c<double>(xpMat, xpMat_G, bin_index, typed_index, typed_bin_index, typed_value, lambda, maf, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
