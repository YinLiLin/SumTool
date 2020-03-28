#include <RcppArmadillo.h>
#include <iostream>
#include <R_ext/Print.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::vec getwt(const double h2, const arma::vec ld, const arma::vec wld, const double N, const int M, const double intercept){
	double h;
	if(h2 <= 0){
		h = 0;
	}else if(h2 >= 1){
		h = 1;
	}else{
		h = h2;
	}
	double coeff = h * N / M;
	arma::vec hetero = 0.5 / square(intercept + coeff * ld);
	arma::vec wt = hetero % (1 / wld);
	return sqrt(wt);
}

arma::mat wtx(const arma::mat x, const arma::vec wt_sqrt){
	arma::vec wt_norm = wt_sqrt / sum(wt_sqrt);
	arma::mat wtx = x.each_col() % wt_norm;
	return wtx;
}

// [[Rcpp::export]]
SEXP ldreg_h2(const arma::vec & z2, const arma::vec & ld, const arma::vec & wld, const int M, const int N, const double maxz2 = 30, const int nblock = 200, const int rep = 3){

	uvec indx = find(z2 < maxz2);
	int m = indx.n_elem;
	if(m < nblock){
		indx = find(z2 <= z2.max());
		m = indx.n_elem;
	}
	double h2 = M * mean(z2 - 1) / (N * mean(ld));
	double h2_se = 0.0;

	// estimating intercept
	arma::mat x(m, 2);
	arma::mat y(m, 1);
	double intercept = 1.0;
	double intercept_se = 0.0;
	x.col(0).fill(intercept);
	x.col(1) = ld(indx);
	y.col(0) = z2(indx);

	arma::mat xwt, ywt, coeff;
	for(int i = 0; i < rep; i++){
		arma::vec wt = getwt(h2, ld(indx), wld(indx), N, M, intercept);
		xwt = wtx(x, wt);
		ywt = wtx(y, wt);
		coeff = inv(xwt.t() * xwt) * (xwt.t() * ywt);
		h2 = coeff(1, 0) * M / N;
		intercept = coeff(0, 0);
	}

	int block_len = m / nblock;
	int start = 0;
	int end = 0;
	arma::mat xtx(2 * nblock, 2);
	arma::mat xty(2, nblock);
	arma::mat xtx_sum(2, 2); xtx_sum.fill(0);
	arma::mat xty_sum(2, 1); xty_sum.fill(0);
	for(int j = 0; j < nblock; j++){
		end = start + block_len - 1;
		if(end >= m || (j == (nblock - 1) && end < (m - 1))){end = m - 1;}
		xtx.rows(j * 2, j * 2 + 1) = xwt.rows(start, end).t() * xwt.rows(start, end);
		xty.col(j) = xwt.rows(start, end).t() * ywt.rows(start, end);
		xty_sum.col(0) += xty.col(j);
		xtx_sum.rows(0, 1) += xtx.rows(j * 2, j * 2 + 1);
		if(end == m - 1){
			break;
		}
		start += block_len;
	}
	coeff = inv(xtx_sum) * xty_sum;

	arma::mat xtx_sum_delete(2, 2);
	arma::mat xty_sum_delete(2, 1);
	arma::mat coeff_jackknife(2, nblock);
	arma::vec h2_jackknife(nblock);
	arma::vec intercept_jackknife(nblock);
	for(int j = 0; j < nblock; j++){
		xtx_sum_delete = xtx_sum - xtx.rows(j * 2, j * 2 + 1);
		xty_sum_delete = xty_sum - xty.col(j);
		coeff_jackknife.col(j) = inv(xtx_sum_delete) * xty_sum_delete;
		intercept_jackknife[j] = nblock * coeff(0, 0) - (nblock - 1) * coeff_jackknife(0, j);
		h2_jackknife[j] = nblock * coeff(1, 0) - (nblock - 1) * coeff_jackknife(1, j);
	}

	intercept = mean(intercept_jackknife);
	intercept_se = stddev(intercept_jackknife) / sqrt(nblock);
	// h2 = mean(h2_jackknife) * M / N;
	// h2_se =  M * stddev(h2_jackknife) / sqrt(nblock) / N;

	// estimating h2
	m = z2.n_elem;
	block_len = m / nblock;
	arma::mat yerr(m, 1);
	yerr.col(0) = z2 - intercept;
	h2 = M * mean(z2 - 1) / (N * mean(ld));

	arma::mat xwt1, ywt1, coeff1;
	for(int i = 0; i < rep; i++){
		arma::vec wt = getwt(h2, ld, wld, N, M, intercept);
		xwt1 = wtx(ld, wt);
		ywt1 = wtx(yerr, wt);
		coeff1 = inv(xwt1.t() * xwt1) * (xwt1.t() * ywt1);
		h2 = coeff1(0, 0) * M / N;
	}

	start = 0;
	end = 0;
	arma::vec xtx1(nblock);
	arma::vec xty1(nblock);
	double xtx_sum1 = 0.0;
	double xty_sum1 = 0.0;
	for(int j = 0; j < nblock; j++){
		end = start + block_len - 1;
		if(end >= m || (j == (nblock - 1) && end < (m - 1))){end = m - 1;}
		arma::mat temp1 = xwt1.rows(start, end).t() * xwt1.rows(start, end);
		xtx1[j] = temp1(0, 0);
		arma::mat temp2 = xwt1.rows(start, end).t() * ywt1.rows(start, end);
		xty1[j] = temp2(0, 0);
		xty_sum1 += xty1[j];
		xtx_sum1 += xtx1[j];
		if(end == m - 1){
			break;
		}
		start += block_len;
	}
	h2 = xty_sum1 / xtx_sum1;

	double xtx_sum_delete1 = 0.0;
	double xty_sum_delete1 = 0.0;
	arma::vec coeff_jackknife1(nblock);
	arma::vec h2_jackknife1(nblock);
	for(int j = 0; j < nblock; j++){
		xtx_sum_delete1 = xtx_sum1 - xtx1[j];
		xty_sum_delete1 = xty_sum1 - xty1[j];
		coeff_jackknife1[j] = xty_sum_delete1 / xtx_sum_delete1;
		h2_jackknife1[j] = nblock * h2 - (nblock - 1) * coeff_jackknife1[j];
	}

	h2 = mean(h2_jackknife1) * M / N;
	h2_se = M * stddev(h2_jackknife1) / sqrt(nblock) / N;

	arma::vec res(4);
	res[2] = intercept;
	res[3] = intercept_se;
	res[0] = h2;
	res[1] = h2_se;

	return wrap(res);
}

// [[Rcpp::export]]
SEXP ldreg_rg(const arma::vec & z_t1, const arma::vec & z_t2, const arma::vec & ld, const arma::vec & wld, const int M, const int N1, const int N2, const int nblock = 200){

	arma::vec zz = z_t1 % z_t2;
	arma::vec zz1 = z_t1 % z_t1;
	arma::vec zz2 = z_t2 % z_t2;
	int m = zz.n_elem;
	double h2 = 0.0, h2_1 = 0.0, h2_2 = 0;
	double h2_se = 0.0;

	// estimating intercept
	arma::mat x(m, 2);
	arma::mat y(m, 1), y1(m, 1), y2(m, 1);
	double intercept = 1.0, intercept1 = 1.0, intercept2 = 1.0;
	double intercept_se = 0.0;
	x.col(0).fill(intercept);
	x.col(1) = ld;
	y.col(0) = zz;
	y1.col(0) = zz1;
	y2.col(0) = zz2;

	arma::mat xwt, ywt, coeff;
	h2_1 = M * mean(zz1 - 1) / (sqrt(N1) * mean(ld));
	h2_2 = M * mean(zz2 - 1) / (sqrt(N2) * mean(ld));

	arma::vec wt1 = getwt(h2_1, ld, wld, N1, M, intercept1);
	// xwt = wtx(x, wt1);
	// ywt = wtx(y1, wt1);
	// coeff = inv(xwt.t() * xwt) * (xwt.t() * ywt);
	// h2_1 = coeff(1, 0) * M / N1;
	// intercept1 = coeff(0, 0);

	arma::vec wt2 = getwt(h2_2, ld, wld, N2, M, intercept2);
	// xwt = wtx(x, wt2);
	// ywt = wtx(y2, wt2);
	// coeff = inv(xwt.t() * xwt) * (xwt.t() * ywt);
	// h2_2 = coeff(1, 0) * M / N2;
	// intercept2 = coeff(0, 0);

	arma::vec wt = wt1 + wt2;
	xwt = wtx(x, wt);
	ywt = wtx(y, wt);
	coeff = inv(xwt.t() * xwt) * (xwt.t() * ywt);
	h2 = coeff(1, 0) * M / sqrt(N1 * N2);
	intercept = coeff(0, 0);

	int block_len = m / nblock;
	int start = 0;
	int end = 0;
	arma::mat xtx(2 * nblock, 2);
	arma::mat xty(2, nblock);
	arma::mat xtx_sum(2, 2); xtx_sum.fill(0);
	arma::mat xty_sum(2, 1); xty_sum.fill(0);
	for(int j = 0; j < nblock; j++){
		end = start + block_len - 1;
		if(end >= m || (j == (nblock - 1) && end < (m - 1))){end = m - 1;}
		xtx.rows(j * 2, j * 2 + 1) = xwt.rows(start, end).t() * xwt.rows(start, end);
		xty.col(j) = xwt.rows(start, end).t() * ywt.rows(start, end);
		xty_sum.col(0) += xty.col(j);
		xtx_sum.rows(0, 1) += xtx.rows(j * 2, j * 2 + 1);
		if(end == m - 1){
			break;
		}
		start += block_len;
	}
	coeff = inv(xtx_sum) * xty_sum;

	arma::mat xtx_sum_delete(2, 2);
	arma::mat xty_sum_delete(2, 1);
	arma::mat coeff_jackknife(2, nblock);
	arma::vec h2_jackknife(nblock);
	arma::vec intercept_jackknife(nblock);
	for(int j = 0; j < nblock; j++){
		xtx_sum_delete = xtx_sum - xtx.rows(j * 2, j * 2 + 1);
		xty_sum_delete = xty_sum - xty.col(j);
		coeff_jackknife.col(j) = inv(xtx_sum_delete) * xty_sum_delete;
		intercept_jackknife[j] = nblock * coeff(0, 0) - (nblock - 1) * coeff_jackknife(0, j);
		h2_jackknife[j] = nblock * coeff(1, 0) - (nblock - 1) * coeff_jackknife(1, j);
	}

	intercept = mean(intercept_jackknife);
	intercept_se = stddev(intercept_jackknife) / sqrt(nblock);
	h2 = mean(h2_jackknife) * M / sqrt(N1 * N2);
	h2_se =  M * stddev(h2_jackknife) / sqrt(nblock) / sqrt(N1 * N2);

	// estimating h2
	// m = zz.n_elem;
	// block_len = m / nblock;
	arma::mat yerr(m, 1);
	yerr.col(0) = zz - intercept;

	arma::mat xwt1 = wtx(ld, wt);
	arma::mat ywt1 = wtx(yerr, wt);
	arma::mat coeff1 = inv(xwt1.t() * xwt1) * (xwt1.t() * ywt1);
	h2 = coeff1(0, 0) * M / sqrt(N1 * N2);

	start = 0;
	end = 0;
	arma::vec xtx1(nblock);
	arma::vec xty1(nblock);
	double xtx_sum1 = 0.0;
	double xty_sum1 = 0.0;
	for(int j = 0; j < nblock; j++){
		end = start + block_len - 1;
		if(end >= m || (j == (nblock - 1) && end < (m - 1))){end = m - 1;}
		arma::mat temp1 = xwt1.rows(start, end).t() * xwt1.rows(start, end);
		xtx1[j] = temp1(0, 0);
		arma::mat temp2 = xwt1.rows(start, end).t() * ywt1.rows(start, end);
		xty1[j] = temp2(0, 0);
		xty_sum1 += xty1[j];
		xtx_sum1 += xtx1[j];
		if(end == m - 1){
			break;
		}
		start += block_len;
	}
	h2 = xty_sum1 / xtx_sum1;

	double xtx_sum_delete1 = 0.0;
	double xty_sum_delete1 = 0.0;
	arma::vec coeff_jackknife1(nblock);
	arma::vec h2_jackknife1(nblock);
	for(int j = 0; j < nblock; j++){
		xtx_sum_delete1 = xtx_sum1 - xtx1[j];
		xty_sum_delete1 = xty_sum1 - xty1[j];
		coeff_jackknife1[j] = xty_sum_delete1 / xtx_sum_delete1;
		h2_jackknife1[j] = nblock * h2 - (nblock - 1) * coeff_jackknife1[j];
	}

	h2 = mean(h2_jackknife1) * M / sqrt(N1 * N2);
	h2_se = M * stddev(h2_jackknife1) / sqrt(nblock) / sqrt(N1 * N2);

	arma::vec res(4);
	res[2] = intercept;
	res[3] = intercept_se;
	res[0] = h2;
	res[1] = h2_se;

	return wrap(res);
}

