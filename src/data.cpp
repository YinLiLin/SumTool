#include <Rcpp.h>
#include "omp_set.h"
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
void rBed_c(std::string bed_file, XPtr<BigMatrix> pMat, const long maxLine, const double NA_C, const bool add = true, const bool impt = true, const bool verbose = true, const int threads = 0){

	string ending = ".bed";
	if (bed_file.length() <= ending.length() ||
		0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending))
		bed_file += ending;

	omp_setup(threads);

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
	code[3] = add ? 0 : 11;
	code[2] = add ? 1 : 12;
	code[1] = static_cast<T>(NA_C);
	code[0] = add ? 2 : 22;

    std::vector<size_t> ggvec(3);
    if(add){
    	ggvec = {0, 1, 2};
    }else{
    	ggvec = {11, 12, 22};
    }

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
	// Progress progress(cond_n, verbose, pb);
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
	// 		progress.increment();
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
		        std::vector<size_t> na_index = {};
		        std::vector<size_t> counts(3);
		        // size_t counts[3] = { 0 };

		        // count allele, record missing index
		       	int max = 0;
		        T major = add ? 0 : 11;

		        if(add){
		        	for (int j = 0; j < nid; j++) {
			            switch(int(mat[i][j])) {
				            case 0: counts[0]++; break;
				            case 1: counts[1]++; break;
				            case 2: counts[2]++; break;
				            default: na_index.push_back(j);
		            	}
		        	}	
		        }else{
		        	for (int j = 0; j < nid; j++) {
			            switch(int(mat[i][j])) {
				            case 11: counts[0]++; break;
				            case 12: counts[1]++; break;
				            case 22: counts[2]++; break;
				            default: na_index.push_back(j);
		            	}
		        	}
		        }

	        	for(size_t j = 0; j < counts.size(); j++){
	        		if(counts[j] > max){
	        			max = counts[j];
	        			major = ggvec[j];
	        		}
	        	}

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
void rBed_c(std::string bfile, const SEXP pBigMat, const long maxLine, const bool add = true, const bool impt = true, const bool verbose = true, const int threads = 0){
	
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return rBed_c<char>(bfile, xpMat, maxLine, NA_CHAR, add, impt, verbose, threads);
	case 2:
		return rBed_c<short>(bfile, xpMat, maxLine, NA_SHORT, add, impt, verbose, threads);
	case 4:
		return rBed_c<int>(bfile, xpMat, maxLine, NA_INTEGER, add, impt, verbose, threads);
	case 8:
		return rBed_c<double>(bfile, xpMat, maxLine, NA_REAL, add, impt, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
void wBed_c(XPtr<BigMatrix> pMat, std::string bed_file, double NA_C, int threads=0, bool verbose=true) {
    
	omp_setup(threads);

    // check input
    string ending = ".bed";
    if (bed_file.length() <= ending.length() ||
        0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending)) {
        bed_file += ending;
    }
    
    // define
    T c;
	int m = pMat->ncol();
	int nid = pMat->nrow();
    int n = pMat->nrow() / 4;  // 4 individual = 1 bit
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
void wBed_c(SEXP pBigMat, std::string bed_file, int threads = 0, bool verbose = true) {
    
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return wBed_c<char>(xpMat, bed_file, NA_CHAR, threads, verbose);
    case 2:
        return wBed_c<short>(xpMat, bed_file, NA_SHORT, threads, verbose);
    case 4:
        return wBed_c<int>(xpMat, bed_file, NA_INTEGER, threads, verbose);
    case 8:
        return wBed_c<double>(xpMat, bed_file, NA_REAL, threads, verbose);
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
//' @export
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
//' @export
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
List VCF_Row_Col(std::string vcf_file){

    string line;
    vector<string> ind;
    int n = 0;
    int m = 0;
    
    ifstream file(vcf_file);
    if (!file) throw Rcpp::exception(("Error: can not open the file [" + vcf_file + "].").c_str());

    // Skip Header
    string prefix("#CHROM");
    bool have_header = false;
    bool record = false;
    while (file) {
        getline(file, line);
        if (!line.compare(0, prefix.size(), prefix)) {
            have_header = true;
            record = true;
            boost::split(ind, line, boost::is_any_of("\t"));
		    vector<string>(ind.begin() + 9, ind.end()).swap(ind);   // delete first 9 columns
		    n = ind.size();
        }
        if(record)	m++;
    }
    if (!have_header) {
        Rcpp::stop("ERROR: Wrong VCF file, no line begin with \"#CHROM\".");
    }
    return List::create(_["n"] = n,
                        _["m"] = m - 2);
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

template <typename T>
List rVCF_c(const std::string vcf_file, XPtr<BigMatrix> pMat, const long maxLine, const std::string out, const double NA_C, const bool impt = true, const bool add = true, const bool verbose = true, const int threads = 0) {

    // Define
    string line;
    vector<string> l;
    int m = pMat->ncol();
	int n = pMat->nrow();
	MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
	vector<string> snp(m);
	vector<string> chr(m);
	vector<string> pos(m);
	vector<string> a1(m);
	vector<string> a2(m);
	NumericVector miss(m); miss.fill(0);
	MinimalProgressBar_plus pb;
	Progress progress(m, verbose, pb);

    omp_setup(threads);

	ifstream file(vcf_file);
    if (!file) throw Rcpp::exception(("Error: can not open the file [" + vcf_file + "].").c_str());

    // Skip Header
    string prefix("#CHROM");
    bool have_header = false;
    while (file) {
        getline(file, line);
        if (!line.compare(0, prefix.size(), prefix)) {
            have_header = true;
            break; 
        }
    }
    if (!have_header) {
        Rcpp::stop("ERROR: Wrong VCF file, no line begin with \"#CHROM\".");
    }

    std::vector<size_t> ggvec(4);
    if(add){
    	ggvec = {0, 1, 1, 2};
    }else{
    	ggvec = {11, 12, 21, 22};
    }
    
    vector<string> buffer;
    int m0 = 0;
    while (file) {
        buffer.clear();
        for (int i = 0; (maxLine <= 0 || i < maxLine); i++) {
        	getline(file, line);
            if (line.length() > 1) {    // Handling the blank line at the end of the file.
                buffer.push_back(line);
            }else{
            	break;
            }
        }

        #pragma omp parallel for private(l)
        for (std::size_t i = 0; i < buffer.size(); i++) {
            
            boost::split(l, buffer[i], boost::is_any_of("\t"));
            
            if (l[2] == ".") {      // snp name missing
           		l[2] = l[0] + '_' + l[1];
        	}
        	snp[m0 + i] = l[2];
        	chr[m0 + i] = l[0];
        	pos[m0 + i] = l[1];
        	a1[m0 + i] = l[3];
        	a2[m0 + i] = l[4];

            vector<string>(l.begin() + 9, l.end()).swap(l);
            for (std::size_t j = 0; j < l.size(); j++) {
            	int gg;
            	if(l[j][0] == '0' && l[j][2] == '0'){
            		gg = ggvec[0];
            	}
            	else if (l[j][0] == '1' && l[j][2] == '1'){
            		gg = ggvec[3];
            	}
            	else if (l[j][0] == '0' && l[j][2] == '1'){
            		gg = ggvec[1];
            	}
            	else if (l[j][0] == '1' && l[j][2] == '0'){
            		gg = ggvec[2];
            	} else {
            		gg = static_cast<T>(NA_C);
					miss[m0 + i] = 1;
            	}
                mat[m0 + i][j] = gg;
            }
            progress.increment();
        }
        m0 += buffer.size();
    }
	file.close();

	ofstream map(out + ".map");
	map << "SNP\tCHROM\tPOS\tREF\tALT" << endl;
	for(int i = 0; i < m; i++){
		map << snp[i] << "\t" << chr[i] << "\t" << pos[i] << "\t" << a1[i] << "\t" << a2[i] << endl;
	}
	map.close();
    
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
		        std::vector<size_t> na_index = {};
		        std::vector<size_t> counts(4);

	        	int max = 0;
	        	T major = add ? 0 : 11;

		        // count allele, record missing index 
		        if(add){
		        	for (int j = 0; j < n; j++) {
			            switch(int(mat[i][j])) {
				            case 0: counts[0]++; break;
				            case 1: counts[1]++; break;
				            case 2: counts[2]++; break;
				            default: na_index.push_back(j);
		            	}
		        	}
		        }else{
		        	for (int j = 0; j < n; j++) {
			            switch(int(mat[i][j])) {
				            case 11: counts[0]++; break;
				            case 12: counts[1]++; break;
				            case 21: counts[2]++; break;
				            case 22: counts[3]++; break;
				            default: na_index.push_back(j);
		            	}
		        	}
		        }

	        	for(size_t j = 0; j < counts.size(); j++){
	        		if(counts[j] > max){
	        			max = counts[j];
	        			major = ggvec[j];
	        		}
	        	}
	        	
	        	// impute
		        for (auto&& x: na_index) {
		            mat[i][x] = major;   
		        }
		    }
		}
	}

	return List::create(
		Named("SNP") = snp,
		Named("Chr") = chr,
		Named("Pos") = pos,
		Named("Ref") = a1,
		Named("Alt") = a2
	);
}

// [[Rcpp::export]]
List rVCF_c(const std::string vcf_file, const SEXP pBigMat, const long maxLine, const std::string out, const bool impt = true, const bool add = true, const bool verbose = true, const int threads = 0){
	
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return rVCF_c<char>(vcf_file, xpMat, maxLine, out, NA_CHAR, impt, add, verbose, threads);
	case 2:
		return rVCF_c<short>(vcf_file, xpMat, maxLine, out, NA_SHORT, impt, add, verbose, threads);
	case 4:
		return rVCF_c<int>(vcf_file, xpMat, maxLine, out, NA_INTEGER, impt, add, verbose, threads);
	case 8:
		return rVCF_c<double>(vcf_file, xpMat, maxLine, out, NA_REAL, impt, add, verbose, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
