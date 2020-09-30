// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// SImputeZ_bin_c
SEXP SImputeZ_bin_c(SEXP pBigMat, const IntegerVector& bin_index, const IntegerVector& typed_index, const arma::mat& typed_value, const double lambda, const double maf, const bool verbose, const int threads);
RcppExport SEXP _SumTool_SImputeZ_bin_c(SEXP pBigMatSEXP, SEXP bin_indexSEXP, SEXP typed_indexSEXP, SEXP typed_valueSEXP, SEXP lambdaSEXP, SEXP mafSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type bin_index(bin_indexSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type typed_index(typed_indexSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type typed_value(typed_valueSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(SImputeZ_bin_c(pBigMat, bin_index, typed_index, typed_value, lambda, maf, verbose, threads));
    return rcpp_result_gen;
END_RCPP
}
// SImputeZ_ld_bin_c
SEXP SImputeZ_ld_bin_c(SEXP pBigMat, SEXP pBigMat_G, const IntegerVector& bin_index, const IntegerVector& typed_index, const IntegerVector& typed_bin_index, const arma::mat& typed_value, const double lambda, const double maf, const bool verbose, const int threads);
RcppExport SEXP _SumTool_SImputeZ_ld_bin_c(SEXP pBigMatSEXP, SEXP pBigMat_GSEXP, SEXP bin_indexSEXP, SEXP typed_indexSEXP, SEXP typed_bin_indexSEXP, SEXP typed_valueSEXP, SEXP lambdaSEXP, SEXP mafSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat_G(pBigMat_GSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type bin_index(bin_indexSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type typed_index(typed_indexSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type typed_bin_index(typed_bin_indexSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type typed_value(typed_valueSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(SImputeZ_ld_bin_c(pBigMat, pBigMat_G, bin_index, typed_index, typed_bin_index, typed_value, lambda, maf, verbose, threads));
    return rcpp_result_gen;
END_RCPP
}
// rBed_c
void rBed_c(std::string bfile, const SEXP pBigMat, const long maxLine, const bool add, const bool impt, const bool verbose, const int threads);
RcppExport SEXP _SumTool_rBed_c(SEXP bfileSEXP, SEXP pBigMatSEXP, SEXP maxLineSEXP, SEXP addSEXP, SEXP imptSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bfile(bfileSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const long >::type maxLine(maxLineSEXP);
    Rcpp::traits::input_parameter< const bool >::type add(addSEXP);
    Rcpp::traits::input_parameter< const bool >::type impt(imptSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rBed_c(bfile, pBigMat, maxLine, add, impt, verbose, threads);
    return R_NilValue;
END_RCPP
}
// wBed_c
void wBed_c(SEXP pBigMat, std::string bed_file, int threads, bool verbose);
RcppExport SEXP _SumTool_wBed_c(SEXP pBigMatSEXP, SEXP bed_fileSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< std::string >::type bed_file(bed_fileSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    wBed_c(pBigMat, bed_file, threads, verbose);
    return R_NilValue;
END_RCPP
}
// FileNcol
int FileNcol(std::string filename);
RcppExport SEXP _SumTool_FileNcol(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(FileNcol(filename));
    return rcpp_result_gen;
END_RCPP
}
// FileNrow
int FileNrow(std::string filename);
RcppExport SEXP _SumTool_FileNrow(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(FileNrow(filename));
    return rcpp_result_gen;
END_RCPP
}
// VCF_Row_Col
List VCF_Row_Col(std::string vcf_file);
RcppExport SEXP _SumTool_VCF_Row_Col(SEXP vcf_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type vcf_file(vcf_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(VCF_Row_Col(vcf_file));
    return rcpp_result_gen;
END_RCPP
}
// rMap_c
List rMap_c(std::string map_file, const Nullable<std::string> out);
RcppExport SEXP _SumTool_rMap_c(SEXP map_fileSEXP, SEXP outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type map_file(map_fileSEXP);
    Rcpp::traits::input_parameter< const Nullable<std::string> >::type out(outSEXP);
    rcpp_result_gen = Rcpp::wrap(rMap_c(map_file, out));
    return rcpp_result_gen;
END_RCPP
}
// rVCF_c
List rVCF_c(const std::string vcf_file, const SEXP pBigMat, const long maxLine, const std::string out, const bool impt, const bool add, const bool verbose, const int threads);
RcppExport SEXP _SumTool_rVCF_c(SEXP vcf_fileSEXP, SEXP pBigMatSEXP, SEXP maxLineSEXP, SEXP outSEXP, SEXP imptSEXP, SEXP addSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type vcf_file(vcf_fileSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const long >::type maxLine(maxLineSEXP);
    Rcpp::traits::input_parameter< const std::string >::type out(outSEXP);
    Rcpp::traits::input_parameter< const bool >::type impt(imptSEXP);
    Rcpp::traits::input_parameter< const bool >::type add(addSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rVCF_c(vcf_file, pBigMat, maxLine, out, impt, add, verbose, threads));
    return rcpp_result_gen;
END_RCPP
}
// SImpute_LD_bigm_c
void SImpute_LD_bigm_c(arma::mat& genomat, SEXP ldbigmat, const Nullable<IntegerVector> index, const int chisq, const double lambda, const bool haps, const int threads, const bool verbose);
RcppExport SEXP _SumTool_SImpute_LD_bigm_c(SEXP genomatSEXP, SEXP ldbigmatSEXP, SEXP indexSEXP, SEXP chisqSEXP, SEXP lambdaSEXP, SEXP hapsSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type genomat(genomatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ldbigmat(ldbigmatSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const int >::type chisq(chisqSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type haps(hapsSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    SImpute_LD_bigm_c(genomat, ldbigmat, index, chisq, lambda, haps, threads, verbose);
    return R_NilValue;
END_RCPP
}
// SImpute_LD_norm_c
SEXP SImpute_LD_norm_c(SEXP pBigMat, const Nullable<IntegerVector> index, const int chisq, const double lambda, const bool haps, const int threads, const bool verbose);
RcppExport SEXP _SumTool_SImpute_LD_norm_c(SEXP pBigMatSEXP, SEXP indexSEXP, SEXP chisqSEXP, SEXP lambdaSEXP, SEXP hapsSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const int >::type chisq(chisqSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type haps(hapsSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(SImpute_LD_norm_c(pBigMat, index, chisq, lambda, haps, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// LDcor_c
arma::mat LDcor_c(SEXP pBigMat1, SEXP pBigMat2, const IntegerVector index1, const IntegerVector index2, const int threads);
RcppExport SEXP _SumTool_LDcor_c(SEXP pBigMat1SEXP, SEXP pBigMat2SEXP, SEXP index1SEXP, SEXP index2SEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat1(pBigMat1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat2(pBigMat2SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index1(index1SEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index2(index2SEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(LDcor_c(pBigMat1, pBigMat2, index1, index2, threads));
    return rcpp_result_gen;
END_RCPP
}
// LDprune_c
arma::uvec LDprune_c(SEXP pBigMat, const IntegerVector index, const double r2_cutoff, const int threads, const bool verbose);
RcppExport SEXP _SumTool_LDprune_c(SEXP pBigMatSEXP, SEXP indexSEXP, SEXP r2_cutoffSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const double >::type r2_cutoff(r2_cutoffSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(LDprune_c(pBigMat, index, r2_cutoff, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// LDclump_c
arma::uvec LDclump_c(SEXP pBigMat, const IntegerVector index, const double r2_cutoff, const arma::vec p, const int threads, const bool verbose);
RcppExport SEXP _SumTool_LDclump_c(SEXP pBigMatSEXP, SEXP indexSEXP, SEXP r2_cutoffSEXP, SEXP pSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const double >::type r2_cutoff(r2_cutoffSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(LDclump_c(pBigMat, index, r2_cutoff, p, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// ldreg_h2
SEXP ldreg_h2(const arma::vec& z2, const arma::vec& ld, const arma::vec& wld, const int M, const int N, const double maxz2, const int nblock, const int rep);
RcppExport SEXP _SumTool_ldreg_h2(SEXP z2SEXP, SEXP ldSEXP, SEXP wldSEXP, SEXP MSEXP, SEXP NSEXP, SEXP maxz2SEXP, SEXP nblockSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type z2(z2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ld(ldSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type wld(wldSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double >::type maxz2(maxz2SEXP);
    Rcpp::traits::input_parameter< const int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< const int >::type rep(repSEXP);
    rcpp_result_gen = Rcpp::wrap(ldreg_h2(z2, ld, wld, M, N, maxz2, nblock, rep));
    return rcpp_result_gen;
END_RCPP
}
// ldreg_rg
SEXP ldreg_rg(const arma::vec& z_t1, const arma::vec& z_t2, const arma::vec& ld, const arma::vec& wld, const int M, const int N1, const int N2, const int nblock);
RcppExport SEXP _SumTool_ldreg_rg(SEXP z_t1SEXP, SEXP z_t2SEXP, SEXP ldSEXP, SEXP wldSEXP, SEXP MSEXP, SEXP N1SEXP, SEXP N2SEXP, SEXP nblockSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type z_t1(z_t1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z_t2(z_t2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ld(ldSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type wld(wldSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type N1(N1SEXP);
    Rcpp::traits::input_parameter< const int >::type N2(N2SEXP);
    Rcpp::traits::input_parameter< const int >::type nblock(nblockSEXP);
    rcpp_result_gen = Rcpp::wrap(ldreg_rg(z_t1, z_t2, ld, wld, M, N1, N2, nblock));
    return rcpp_result_gen;
END_RCPP
}
// LDscore_c
SEXP LDscore_c(SEXP pBigMat, const IntegerVector index, const bool r2, const bool adjust, const int threads, const bool verbose);
RcppExport SEXP _SumTool_LDscore_c(SEXP pBigMatSEXP, SEXP indexSEXP, SEXP r2SEXP, SEXP adjustSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const bool >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< const bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(LDscore_c(pBigMat, index, r2, adjust, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// getCol
NumericVector getCol(SEXP pBigMat, const int col);
RcppExport SEXP _SumTool_getCol(SEXP pBigMatSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const int >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(getCol(pBigMat, col));
    return rcpp_result_gen;
END_RCPP
}
// glm_c
SEXP glm_c(const arma::vec y, const arma::mat X, const IntegerVector indx, SEXP pBigMat, const bool verbose, const int threads);
RcppExport SEXP _SumTool_glm_c(SEXP ySEXP, SEXP XSEXP, SEXP indxSEXP, SEXP pBigMatSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type indx(indxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(glm_c(y, X, indx, pBigMat, verbose, threads));
    return rcpp_result_gen;
END_RCPP
}
// sblup_bin
SEXP sblup_bin(SEXP pBigMat, const double n_gwas, const IntegerVector typed_index, const arma::vec typed_value, const double lambda, const bool verbose, const int threads);
RcppExport SEXP _SumTool_sblup_bin(SEXP pBigMatSEXP, SEXP n_gwasSEXP, SEXP typed_indexSEXP, SEXP typed_valueSEXP, SEXP lambdaSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const double >::type n_gwas(n_gwasSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type typed_index(typed_indexSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type typed_value(typed_valueSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sblup_bin(pBigMat, n_gwas, typed_index, typed_value, lambda, verbose, threads));
    return rcpp_result_gen;
END_RCPP
}
// hasNA
bool hasNA(SEXP pBigMat, const int threads);
RcppExport SEXP _SumTool_hasNA(SEXP pBigMatSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(hasNA(pBigMat, threads));
    return rcpp_result_gen;
END_RCPP
}
// freq_snp
NumericVector freq_snp(SEXP pBigMat, const IntegerVector index_, const int threads);
RcppExport SEXP _SumTool_freq_snp(SEXP pBigMatSEXP, SEXP index_SEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index_(index_SEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(freq_snp(pBigMat, index_, threads));
    return rcpp_result_gen;
END_RCPP
}
// freq_hap
double freq_hap(SEXP pBigMat, const int indx_1, const int indx_2);
RcppExport SEXP _SumTool_freq_hap(SEXP pBigMatSEXP, SEXP indx_1SEXP, SEXP indx_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const int >::type indx_1(indx_1SEXP);
    Rcpp::traits::input_parameter< const int >::type indx_2(indx_2SEXP);
    rcpp_result_gen = Rcpp::wrap(freq_hap(pBigMat, indx_1, indx_2));
    return rcpp_result_gen;
END_RCPP
}
// freq_s
NumericVector freq_s(SEXP pBigMat, const IntegerVector index_, const int threads);
RcppExport SEXP _SumTool_freq_s(SEXP pBigMatSEXP, SEXP index_SEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index_(index_SEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(freq_s(pBigMat, index_, threads));
    return rcpp_result_gen;
END_RCPP
}
// freq_h
double freq_h(SEXP pBigMat, const int indx_1, const int indx_2);
RcppExport SEXP _SumTool_freq_h(SEXP pBigMatSEXP, SEXP indx_1SEXP, SEXP indx_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const int >::type indx_1(indx_1SEXP);
    Rcpp::traits::input_parameter< const int >::type indx_2(indx_2SEXP);
    rcpp_result_gen = Rcpp::wrap(freq_h(pBigMat, indx_1, indx_2));
    return rcpp_result_gen;
END_RCPP
}
// BigStat
SEXP BigStat(SEXP pBigMat, const IntegerVector index_, const int threads);
RcppExport SEXP _SumTool_BigStat(SEXP pBigMatSEXP, SEXP index_SEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type index_(index_SEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(BigStat(pBigMat, index_, threads));
    return rcpp_result_gen;
END_RCPP
}
// which_c
IntegerVector which_c(const NumericVector x, const double value, const int c);
RcppExport SEXP _SumTool_which_c(SEXP xSEXP, SEXP valueSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type value(valueSEXP);
    Rcpp::traits::input_parameter< const int >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(which_c(x, value, c));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SumTool_SImputeZ_bin_c", (DL_FUNC) &_SumTool_SImputeZ_bin_c, 8},
    {"_SumTool_SImputeZ_ld_bin_c", (DL_FUNC) &_SumTool_SImputeZ_ld_bin_c, 10},
    {"_SumTool_rBed_c", (DL_FUNC) &_SumTool_rBed_c, 7},
    {"_SumTool_wBed_c", (DL_FUNC) &_SumTool_wBed_c, 4},
    {"_SumTool_FileNcol", (DL_FUNC) &_SumTool_FileNcol, 1},
    {"_SumTool_FileNrow", (DL_FUNC) &_SumTool_FileNrow, 1},
    {"_SumTool_VCF_Row_Col", (DL_FUNC) &_SumTool_VCF_Row_Col, 1},
    {"_SumTool_rMap_c", (DL_FUNC) &_SumTool_rMap_c, 2},
    {"_SumTool_rVCF_c", (DL_FUNC) &_SumTool_rVCF_c, 8},
    {"_SumTool_SImpute_LD_bigm_c", (DL_FUNC) &_SumTool_SImpute_LD_bigm_c, 8},
    {"_SumTool_SImpute_LD_norm_c", (DL_FUNC) &_SumTool_SImpute_LD_norm_c, 7},
    {"_SumTool_LDcor_c", (DL_FUNC) &_SumTool_LDcor_c, 5},
    {"_SumTool_LDprune_c", (DL_FUNC) &_SumTool_LDprune_c, 5},
    {"_SumTool_LDclump_c", (DL_FUNC) &_SumTool_LDclump_c, 6},
    {"_SumTool_ldreg_h2", (DL_FUNC) &_SumTool_ldreg_h2, 8},
    {"_SumTool_ldreg_rg", (DL_FUNC) &_SumTool_ldreg_rg, 8},
    {"_SumTool_LDscore_c", (DL_FUNC) &_SumTool_LDscore_c, 6},
    {"_SumTool_getCol", (DL_FUNC) &_SumTool_getCol, 2},
    {"_SumTool_glm_c", (DL_FUNC) &_SumTool_glm_c, 6},
    {"_SumTool_sblup_bin", (DL_FUNC) &_SumTool_sblup_bin, 7},
    {"_SumTool_hasNA", (DL_FUNC) &_SumTool_hasNA, 2},
    {"_SumTool_freq_snp", (DL_FUNC) &_SumTool_freq_snp, 3},
    {"_SumTool_freq_hap", (DL_FUNC) &_SumTool_freq_hap, 3},
    {"_SumTool_freq_s", (DL_FUNC) &_SumTool_freq_s, 3},
    {"_SumTool_freq_h", (DL_FUNC) &_SumTool_freq_h, 3},
    {"_SumTool_BigStat", (DL_FUNC) &_SumTool_BigStat, 3},
    {"_SumTool_which_c", (DL_FUNC) &_SumTool_which_c, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SumTool(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
