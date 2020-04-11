#include <Rcpp.h>
#include <omp.h>
#include <fstream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include "probar.h"
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(bigmemory, BH)]]

using namespace std;
using namespace Rcpp;

template <typename T>
void rData_c(std::string bed_file, XPtr<BigMatrix> pMat, const long maxLine, const double NA_C, const bool impt = true, const bool verbose = true, const int threads = 0){

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
	NumericVector miss(m); miss.fill(0);
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
	if (!fin) throw Rcpp::exception(("Error: can not open the file [" + bed_file + "].").c_str());
	fseek(fin, 0, SEEK_END);
	long length = ftell(fin);
	rewind(fin);

	// get buffer_size
	buffer_size = maxLine > 0 ? (maxLine * n) : (length - 3);

	int n_block = (length - 3) / buffer_size;
	if ((length - 3) % buffer_size != 0) { n_block++; }

	buffer = new char [3];
	fread(buffer, 1, 3, fin);

	int r, c, x;
	uint8_t p;

// 	int cond_n = 0;
// 	for (int i = 0; i < n_block; i++) {
// 		cond_n += min(buffer_size, (length - 3 - i * buffer_size));
// 	}

	MinimalProgressBar_plus pb;
	// Progress progress(cond_n, verbose);
// 	Progress progress(cond_n, verbose, pb);
	// Rcout << sum(cond) << endl;
	Progress progress(n_block, verbose, pb);

	long block_start;
	for (int i = 0; i < n_block; i++) {
		buffer = new char [buffer_size];
		fread(buffer, 1, buffer_size, fin);

		block_start = i * buffer_size;
		int cond = min(buffer_size, (length - 3 - block_start));

		#pragma omp parallel for schedule(dynamic) private(r, c, p, x)
		for (int j = 0; j < cond; j++) {
			// bit -> item in matrix
			r = (block_start + j) / n;
			c = (block_start + j) % n * 4;
			p = buffer[j];
			for (x = 0; x < 4 && (c + x) < nid; x++) {
				T gg = code[(p >> (2*x)) & 0x03];
				mat[r][c + x] = gg;
				if(gg == NA_C){
					miss[r] = 1;
				}
			}
// 			progress.increment();
		}
		progress.increment();
	}
	fclose(fin);

	int NMISS = 0;
	for (int i = 0; i < m; i++) {
		if(miss[i]){
			Rcout << "Warning: Missing values in genotype exist!" << endl;
			break;
		}else{
			NMISS++;
		}
	}

	if(impt && (NMISS != m)){

		Rcout << "Imputing missing values by major genotype..." << endl;

	 	// impute
	 	#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < m; i++) {
			if(miss[i]){
		        std::vector<size_t> na_index = {};;
		        size_t counts[3] = { 0 };

		        // count allele, record missing index 
	        	for (int j = 0; j < nid; j++) {
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
	}

	return;
}

// [[Rcpp::export]]
void rData_c(std::string bfile, const SEXP pBigMat, const long maxLine, const bool impt = true, const bool verbose = true, const int threads = 0){
	
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return rData_c<char>(bfile, xpMat, maxLine, NA_CHAR, impt, verbose, threads);
	case 2:
		return rData_c<short>(bfile, xpMat, maxLine, NA_SHORT, impt, verbose, threads);
	case 4:
		return rData_c<int>(bfile, xpMat, maxLine, NA_INTEGER, impt, verbose, threads);
	case 8:
		return rData_c<double>(bfile, xpMat, maxLine, NA_REAL, impt, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
void wData_c(XPtr<BigMatrix> pMat, std::string bed_file, double NA_C, int threads=0, bool verbose=true) {
    
	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

    // check input
    string ending = ".bed";
    if (bed_file.length() <= ending.length() ||
        0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending)) {
        bed_file += ending;
    }
    
    // define
    T c;
	long m = pMat->ncol();
	long nid = pMat->nrow();
    long n = pMat->nrow() / 4;  // 4 individual = 1 bit
    if (n % 4 != 0) 
        n++;
    
    vector<uint8_t> geno(n);
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    FILE *fout;
    fout = fopen(bed_file.c_str(), "wb");
    
    // progress bar
   	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);

    // Progress progress(m, verbose);
    
    // magic number of bfile
    const unsigned char magic_bytes[] = { 0x6c, 0x1b, 0x01 };
    fwrite((char*)magic_bytes, 1, 3, fout);
    
    // map
    std::map<T, int> code;
    code[0] = 3;
    code[1] = 2;
    code[2] = 0;
    code[static_cast<T>(NA_C)] = 1;
    
    // write bfile
    for (int i = 0; i < m; i++) {
        #pragma omp parallel for schedule(dynamic) private(c)
        for (int j = 0; j < n; j++) {
            uint8_t p = 0;
            for (int x = 0; x < 4 && (4 * j + x) < nid; x++) {
                c = mat[i][(4 * j + x)];
                p |= code[c] << (x*2);
            }
            geno[j] = p;
        }
        fwrite((char*)geno.data(), 1, geno.size(), fout);
        progress.increment();
    }
    fclose(fout);
    return;
}

// [[Rcpp::export]]
void wData_c(SEXP pBigMat, std::string bed_file, int threads = 0, bool verbose = true) {
    
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return wData_c<char>(xpMat, bed_file, NA_CHAR, threads, verbose);
    case 2:
        return wData_c<short>(xpMat, bed_file, NA_SHORT, threads, verbose);
    case 4:
        return wData_c<int>(xpMat, bed_file, NA_INTEGER, threads, verbose);
    case 8:
        return wData_c<double>(xpMat, bed_file, NA_REAL, threads, verbose);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

//' Column number count
//'
//' To count number of column of a file
//'
//' @param filename, name of a file.
//' @examples
//' # reading data
//' bim_path <- system.file("extdata", "ref_geno.bim", package = "SumTool")
//' n <- FileNcol(bim_path)
// [[Rcpp::export]]
int FileNcol(std::string filename) {
    // Define
    string line;
    vector<string> l;
    ifstream file(filename);
    
    if (!file) throw Rcpp::exception(("Error: can not open the file [" + filename + "].").c_str());

    getline(file, line);
    boost::split(l, line, boost::is_any_of("\t ,"));
    int m = l.size();

    file.close();
    return m;
}

//' Row number count
//'
//' To count number of row of a file
//'
//' @param filename, name of a file.
//' @examples
//' # reading data
//' bim_path <- system.file("extdata", "ref_geno.bim", package = "SumTool")
//' n <- FileNrow(bim_path)
// [[Rcpp::export]]
int FileNrow(std::string filename) {
    // Define
    string line;
    int m = 0; 
    ifstream file(filename);
    
    if (!file) throw Rcpp::exception(("Error: can not open the file [" + filename + "].").c_str());

    while (getline(file, line))
        m++;
    
    file.close();
    return m;
}

// [[Rcpp::export]]
List rMap_c(std::string map_file, const Nullable<std::string> out = R_NilValue){

	bool fileout;
	std::string filename;
	if(out.isNotNull()){
		fileout = true;
		filename = as<std::string>(out);
	}else{
		fileout = false;
	}

	int n = FileNrow(map_file);
	string line;
	vector<string> l;
	vector<string> snp(n);
	vector<string> chr(n);
	vector<string> pos(n);
	vector<string> a1(n);
	vector<string> a2(n);

	ifstream file(map_file);
	if (!file) throw Rcpp::exception(("Error: can not open the file [" + map_file + "].").c_str());

	int idx = 0;
	string s, r, a, c, g, p;
	// float g, p;
	while (file >> c >> s >> g >> p >> r >> a) {

		snp[idx] = s;
		chr[idx] = c;
		pos[idx] = p;
		a1[idx] = r;
		a2[idx] = a;

		idx++;
	}
	// while (getline(file, line)) {

	// 	boost::split(l, line, boost::is_any_of("\t ,"));
	// 	// boost::split(l, line, "\t");

	// 	snp[idx] = l[1];
	// 	chr[idx] = l[0];
	// 	pos[idx] = l[3];
	// 	a1[idx] = l[4];
	// 	a2[idx] = l[5];

	// 	idx++;
	// }
	file.close();

	if(fileout){
		ofstream map(filename + ".map");
		map << "SNP\tCHROM\tPOS\tREF\tALT" << endl;
		for(int i = 0; i < n; i++){
			map << snp[i] << "\t" << chr[i] << "\t" << pos[i] << "\t" << a1[i] << "\t" << a2[i] << endl;
		}
		map.close();
	}

	return List::create( Named("SNP") = snp,
                       Named("Chr") = chr,
                       Named("Pos") = pos,
                       Named("Ref") = a1,
                       Named("Alt") = a2
	);
}
