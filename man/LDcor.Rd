% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SumTool.r
\name{LDcor}
\alias{LDcor}
\title{LD correlation}
\usage{
LDcor(
  geno1 = NULL,
  map1 = NULL,
  geno2 = NULL,
  map2 = NULL,
  w = 1e+06,
  threads = 1,
  verbose = TRUE
)
}
\arguments{
\item{geno1}{big.matrix (n1 * m1), genotype}

\item{map1}{data.frame, the genomic information of SNPs. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2")}

\item{geno2}{big.matrix (n2 * m2), genotype}

\item{map2}{data.frame, the genomic information of SNPs. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2")}

\item{w}{int, size of windows in bp. Default is 1e6}

\item{threads}{number, the number of used threads for parallel process}

\item{verbose}{logical, whether to print the log information}
}
\description{
To compare the LD structure between two populations
}
\examples{
ref_file_path <- system.file("extdata", "ref_geno", package = "SumTool")
gwas_file_path <- system.file("extdata", "gwas_geno", package = "SumTool")
data1 <- read_binary(bfile=ref_file_path, threads=1, verbose=FALSE, out=tempfile())
ref.geno <- data1$geno
ref.map <- data1$map
data2 <- read_binary(bfile=gwas_file_path, threads=1, verbose=FALSE, out=tempfile())
gwas.geno <- data2$geno
gwas.map <- data2$map
ldcor <- LDcor(geno1=ref.geno, map1=ref.map, geno2=gwas.geno, map2=gwas.map, threads=10, verbose=TRUE)
}
