#include <Rcpp.h>
#include <omp.h>
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <progress.hpp>
#include "progress_bar.hpp"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;

template <typename T>
void Data_c(std::string bed_file, XPtr<BigMatrix> pMat, const long maxLine, const double NA_C, const bool verbose = true, const int threads = 0){

	string ending = ".bed";
	if (bed_file.length() <= ending.length() ||
		0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending))
		bed_file += ending;

	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

	long n = pMat->nrow() / 4;  // 4 individual = 1 bit
	int m = pMat->ncol();
	int nid = pMat->nrow();
	NumericVector miss(m); 
	if (pMat->nrow() % 4 != 0)
		n++;
	char *buffer;
	long buffer_size;
	MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);

	// map
	std::map<int, T> code;
	code[3] = 0;
	code[2] = 1;
	code[1] = static_cast<T>(NA_C);
	code[0] = 2;

	// open file
	FILE *fin;
	fin = fopen(bed_file.c_str(), "rb");
	fseek(fin, 0, SEEK_END);
	long length = ftell(fin);
	rewind(fin);

	// get buffer_size
	buffer_size = maxLine > 0 ? (maxLine * n) : (length - 3);

	int n_block = (length - 3) / buffer_size;
	if ((length - 3) % buffer_size != 0) { n_block++; }
	Progress progress(n_block, verbose);

	buffer = new char [3];
	fread(buffer, 1, 3, fin);

	size_t r, c, x;
	uint8_t p;

	long block_start;
	for (int i = 0; i < n_block; i++) {
		buffer = new char [buffer_size];
		fread(buffer, 1, buffer_size, fin);

		block_start = i * buffer_size;
		int cond = min(buffer_size, (length - 3 - i * buffer_size));

		#pragma omp parallel for schedule(dynamic) private(r, c, p, x)
		for (size_t j = 0; j < cond; j++) {
			// bit -> item in matrix
			r = (block_start + j) / n;
			c = (block_start + j) % n * 4;
			p = buffer[j];
			for (x = 0; x < 4 && (c + x) < pMat->nrow(); x++) {
				mat[r][c + x] = code[(p >> (2*x)) & 0x03];
				if(mat[r][c + x] == NA_C){
					miss[r] == 1;
				}
			}
		}
		progress.increment();
	}
	fclose(fin);

	// impute
	for (size_t i = 0; i < m; i++) {
		if(miss[i]){
			Rcout << "(Missings exist! Would be imputed by the major genotype of each SNPs!)" << endl;
			break;
		}
	}

 	#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < m; i++) {
		if(miss[i]){
	        std::vector<size_t> na_index = {};;
	        size_t counts[3] = { 0 };

	        // count allele, record missing index 
        	for (size_t j = 0; j < nid; j++) {
	            switch(int(mat[i][j])) {
		            case 0: counts[0]++; break;
		            case 1: counts[1]++; break;
		            case 2: counts[2]++; break;
		            default: na_index.push_back(j);
            	}
        	}

        	// find major allele
        	T major = counts[2] > counts[1] ? (counts[2] > counts[0] ? 2 : 0) : (counts[1] > counts[0] ? 1 : 0);	
        	
        	// impute
	        for (auto&& x: na_index) {
	            mat[i][x] = major;   
	        }
	    }
	}
	return;
}

// [[Rcpp::export]]
void Data_c(std::string bfile, const SEXP pBigMat, const long maxLine, const bool verbose = true, const int threads = 0){
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return Data_c<char>(bfile, xpMat, maxLine, NA_CHAR, verbose, threads);
	case 2:
		return Data_c<short>(bfile, xpMat, maxLine, NA_SHORT, verbose, threads);
	case 4:
		return Data_c<int>(bfile, xpMat, maxLine, NA_INTEGER, verbose, threads);
	case 8:
		return Data_c<double>(bfile, xpMat, maxLine, NA_REAL, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
