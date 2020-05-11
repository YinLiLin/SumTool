#include "stats.h"
#include "probar.h"

arma::mat GInv(const arma::mat A){
	
	arma::mat ginv;
	if(A.n_rows == 1){
		ginv = 1 / A;
	}else{
		arma::mat U;
		arma::vec s;
		arma::mat V;
		double tol = sqrt(datum::eps);
		
		svd(U,s,V,A);
		U = conv_to<mat>::from(conj(conv_to<cx_mat>::from(U)));
		arma::vec sMax(2); sMax.fill(0);
		sMax[1] = tol * s[0];
		arma::uvec Positive = find(s > sMax.max());
		arma::mat Up = U.cols(Positive);
		Up.each_row() %= 1/s(Positive).t();
		ginv = V.cols(Positive) * Up.t();
	}
	return ginv;
}

template <typename T>
NumericVector getCol(XPtr<BigMatrix> pMat, const int col){
	
	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();

	NumericVector snp(ind);

	for(int i = 0; i < ind; i++){
		snp[i] = genomat[col][i];
	}

	return snp;
}

// [[Rcpp::export]]
NumericVector getCol(SEXP pBigMat, const int col){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return getCol<char>(xpMat, col);
	case 2:
		return getCol<short>(xpMat, col);
	case 4:
		return getCol<int>(xpMat, col);
	case 8:
		return getCol<double>(xpMat, col);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP glm_c(const arma::vec y, const arma::mat X, const IntegerVector indx, XPtr<BigMatrix> pMat, const bool verbose = true, const int threads = 0){
	
	omp_setup(threads);
	
	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int ind = indx.size();
	int mkr = pMat->ncol();
	int q0 = X.n_cols;

	int y_len = y.n_elem;
	if(y_len != ind)
		throw Rcpp::exception("number of individuals not match.!");

	MinimalProgressBar_plus pb;
	Progress progress(mkr, verbose, pb);

	arma::mat iXX = GInv(X.t() * X);
	arma::mat xy = X.t() * y;
	double yy = dot(y, y);
	double maf = 0.0;
	arma::mat res(mkr, 4);
	arma::vec snp(ind);
	arma::mat iXXs(q0 + 1, q0 + 1);

	#pragma omp parallel for schedule(dynamic) private(maf) firstprivate(snp, iXXs)
	for(int i = 0; i < mkr; i++){

		for(int ii = 0; ii < ind; ii++){
			snp[ii] = genomat[i][indx[ii]];
		}
		
		maf = sum(snp) / (2 * ind);
		if(maf > 0.5)	maf = 1 - maf;
		double sy = dot(snp, y);
		double ss = dot(snp, snp);
		arma::mat xs = X.t() * snp;
		arma::mat B21 = xs.t() * iXX;
		double t2 = as_scalar(B21 * xs);
		double invB22 = 1 / (ss - t2);
		arma::mat NeginvB22B21 = -1 * invB22 * B21;

		iXXs(q0, q0)=invB22;

		iXXs.submat(0, 0, q0 - 1, q0 - 1) = iXX + invB22 * B21.t() * B21;
		iXXs(q0, span(0, q0 - 1)) = NeginvB22B21;
		iXXs(span(0, q0 - 1), q0) = NeginvB22B21.t();

        // statistics
        arma::mat rhs(xy.n_rows + 1, 1);
        rhs.rows(0, xy.n_rows - 1) = xy;
        rhs(xy.n_rows, 0) = sy;
		arma::mat beta = iXXs * rhs;
        int df = ind - q0 - 1;
        double ve = (yy - as_scalar(beta.t() * rhs)) / df;

        res(i, 0) = maf;
        res(i, 1) = beta(q0, 0);
        res(i, 2) = sqrt(iXXs(q0, q0) * ve); 
        res(i, 3) = 2 * R::pt(abs(res(i, 1) / res(i, 2)), df, false, false);

        progress.increment();
	}

	return wrap(res);
}

// [[Rcpp::export]]
SEXP glm_c(const arma::vec y, const arma::mat X, const IntegerVector indx, SEXP pBigMat, const bool verbose = true, const int threads = 0){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return glm_c<char>(y, X, indx, xpMat, verbose, threads);
	case 2:
		return glm_c<short>(y, X, indx, xpMat, verbose, threads);
	case 4:
		return glm_c<int>(y, X, indx, xpMat, verbose, threads);
	case 8:
		return glm_c<double>(y, X, indx, xpMat, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
