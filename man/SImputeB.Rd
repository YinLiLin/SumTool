% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SumTool.r
\name{SImputeB}
\alias{SImputeB}
\title{Summary data based imputation on marginal effect}
\usage{
SImputeB(
  ref.geno = NULL,
  ref.map = NULL,
  typed.geno = NULL,
  typed = NULL,
  w = 1e+06,
  b = 5e+05,
  lambda = NULL,
  maf = 1e-06,
  verbose = TRUE,
  threads = 1
)
}
\arguments{
\item{ref.geno}{big.matrix (n1 * m1), reference genotype panel, can be loaded from both binary or vcf file, we recommend using phased vcf file for summary data imputation}

\item{ref.map}{matrix (m1 * 5): SNPs, Chr, position, A1, A2}

\item{typed.geno}{big.matrix (n2 * m2), individual level genotype for typed SNPs (this file is optional)}

\item{typed}{matrix (m2 * 8): SNPs, Chr, position, A1, A2, BETA, SE, N}

\item{w}{number, set the window size in bp, default 1000000}

\item{b}{number, set the buffer size in bp of each window, default 250000}

\item{lambda}{number, ridge regression value on LD matrix of typed SNPs: solve(Rtt + diag(lambda))}

\item{maf}{number, SNPs whose minor allele frequency are lower than set value will not be imputed}

\item{verbose}{logical, whether to print the log information}

\item{threads}{number, the number of used threads for parallel process}
}
\description{
To impute marginal effect using summary data by SImpute/SImpute-LD
}
\examples{
#----------------SImpute----------------#

# get path of the attached files in package
ref_file_path <- system.file("extdata", "ref_geno", package = "SumTool")
typed_b_path <- system.file("extdata", "typed.beta", package = "SumTool")

# reading data
data <- read_binary(bfile=ref_file_path, additive=FALSE, out=tempfile(), 
	threads=1, verbose=FALSE)
ref.geno <- data$geno
ref.map <- data$map
typed_b <- read.table(typed_b_path, header=TRUE)

# Impute marginal effect and se
xx <- SImputeB(ref.geno=ref.geno, ref.map=ref.map, typed=typed_b, threads=1, verbose=FALSE)

#--------------SImpute-LD---------------#

gwas_file_path <- system.file("extdata", "gwas_geno", package = "SumTool")
gwas <- read_binary(bfile=gwas_file_path,  additive=FALSE, out=tempfile(), 
	threads=1, verbose=FALSE)
typed.geno <- gwas$geno
typed.map <- gwas$map
#NOTE: the order of SNPs in 'typed.geno' should be consistent with the order in 'typed_b'.
typed.geno <- deepcopy(typed.geno, cols = match(typed_b[, 1], typed.map[, 1]))

# Impute marginal effect and se
xx <- SImputeB(ref.geno=ref.geno, ref.map=ref.map, typed=typed_b, 
typed.geno=typed.geno, threads=1, verbose=FALSE)
}
