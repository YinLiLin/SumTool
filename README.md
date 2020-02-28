# SumTool [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/YinLiLin/SumTool/issues) [![](https://img.shields.io/badge/Release-v0.99.5-darkred.svg)](https://github.com/YinLiLin/SumTool)

## *A memory-efficient, parallel-accelerated tool for Summary data based Analysis*

Overview
-----
```SumTool``` was designed to implement GWAS summary statistics analysis for big data. It can do ld, ldscore computation, Zscore, Marginal effect imputation, summary data based SNP BLUP (SBLUP) solution, as well as LD regression for heritability and genetic correlation estimation. It featured with memory efficiency and parallel calculation. By the aid of package '[bigmemory](https://cran.r-project.org/web/packages/bigmemory/bigmemory.pdf)', SumTool constructs memory-mapped files for genotype panel on disk instead of loading it all into Random Access Memory (RAM), which makes it possible to handle very big data with limited computation resources. Additionally, all parallel processes are accelerated as fast as possible by OpenMP technology. The functions of SumTool will keep on being enriched with more features based on users feedbacks.

```SumTool``` was developed by [Lilin Yin](https://github.com/YinLiLin) with the support of [Jian Zeng](https://scholar.google.com/citations?user=mOyykToAAAAJ&hl=en) and [Jian Yang](https://scholar.google.com.au/citations?user=aLuqQs8AAAAJ&hl=en). If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/YinLiLin/SumTool/issues/new):point_left:.

Features
-----
- LD 
  - [Linkage disequilibrium matrix](#linkage-disequilibrium)
  - [Linkage disequilibrium score](#ld-score)
  - [Linkage disequilibrium pruning](#ld-pruning)/[clumping](#ld-clumping)
- Impute
  - [Zscore imputation](#impute-zscore)
  - [Beta and SE imputation](#impute-marginal-effect)
- LD regression
  - [Heritability](#estimate-h2)
  - [Genetic correlation](#estimate-rg)
- SBLUP
  - [Joint effect](#estimate-joint-effect)
--- 

Installation
-----
Please install ```devtools``` prior to installing ```SumTool```:
```r
> devtools::install_github("YinLiLin/SumTool")
```
After installed successfully, type ```library(SumTool)``` to use, then type ```lsf.str("package:SumTool")``` see details of all available functions and parameters.

Genotype Converting
-----
With the developing of advaced sequence technologies, the unprecedentedly increased markers across genome make a big challenge in relevant genetic analysis. Loading the genotype into memory directly is highly limited by the computation resources, which becomes the bottleneck of the most of softwares or pipelines. Here we provided a function ```read_plink``` developed by aid of bigmemory package to construct memory-mapped files ('big.matrix') on disk from plink binary files. Instead of reading all genotype data into RAM, it greatly reduce memory cost without significantly increase of computation time.
```r
ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
data <- read_plink(bfile=ref_bfile_path, threads=1)
# bfile: the prefix of binary files
ref.geno <- data$geno
ref.map <- data$map
```
By default, the memory-mapped files are directed into R tempary folder, users could redirect to new path as following:
```r
data <- read_plink(bfile=ref_bfile_path, backingpath="./", descriptorfile="test.desc", backingfile="test.bin", threads=1)

# directly use for the next time
ref.geno <- attach.big.matrix("./test.desc")
```
Linkage Disequilibrium
-----
Here, we simply calculate LD by the Pearson correlation (r) of pairs of SNPs
```r
ld <- LDcal(geno=geno, threads=1)
```
By default, it returns a standard square R matrix with dimension of m by m (m is the number of SNPs), which may cost huge space in memory with the increasing SNPs, users can direct a file on disk by parameter 'out' to store the LD matrix in 'big.matrix' format.

LD Score
-----
LD score is defined as the sum of LD r2 between a variant and all the variants in a region. 
```r
ldscore <- LDsore(geno = geno, map = map, w = 100000, b=50000, threads = 1)
```
In ```LDscore```, we provide a parameter 'r2', users could determine to calculate r or r2. By default, the 'r2' is adjusted by r2adj = r2 - [(1 - r2) / (n -2)], as well as 'r'.

LD Pruning
-----

LD Clumping
-----

Impute Zscore
-----
For Zscore imputation, at least 6 columns should be provided in same order with the example for typed SNPs. No need to separate genome into chromosomes
```r
# get the path of attached example data on SumTool
ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
typed_z_path <- system.file("extdata", "typed.zscore", package = "SumTool")

# reading data
data <- read_plink(bfile=ref_bfile_path, threads=1)
ref.geno <- data$geno
ref.map <- data$map
typed_z <- read.table(typed_z_path, header=TRUE)
head(typed_z)     
         SNP Chr     BP A1 A2     Zscore
1 rs12564807   1 734462  G  A  1.7849800
2  rs3094315   1 752566  G  A  1.3555500
3  rs3131972   1 752721  A  G  1.0441300
4  rs3131969   1 754182  A  G -0.0894082
5  rs1048488   1 760912  C  T -0.1412480

# Impute Zscore
xx <- SImputeZ(ref.geno=ref.geno, ref.map=ref.map, typed=typed_z, w=1000000, threads=1)
```
***NOTE***: For multiple traits, Zscore could be listed in the following columns respectively.<br>
If the individual genotype of summary statistics is available, the imputation accuracy can be improved by using the LD matrix derived from individual genotype rather than reference panel for typed SNPs. 
```r
gwas_bfile_path <- system.file("extdata", "gwas_geno", package = "SumTool")
gwas <- read_plink(bfile=gwas_bfile_path, threads=1)
typed.geno <- gwas$geno
typed.map <- gwas$map
# NOTE: the order of SNPs in 'typed.geno' should be consistent with the order in 'typed_z'.
typed.geno <- deepcopy(typed.geno, cols = match(typed_z[, 1], typed.map[, 1]))
xx <- SImputeZ(ref.geno=ref.geno, ref.map=ref.map, typed=typed_z, typed.geno=typed.geno, w=1000000, threads=1)
```

Impute Marginal Effect
-----
For BETA and SE imputation, limited 8 columns should be provided in same order with the example for typed SNPs. No need to separate genome into chromosomes.<br>
***NOTE***: It is not supported to impute multiple traits at a time.
```r
# get the path of attached example data on SumTool
ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
typed_beta_path <- system.file("extdata", "typed.marginal", package = "SumTool")

# reading data
data <- read_plink(bfile=ref_bfile_path, threads=1)
ref.geno <- data$geno
ref.map <- data$map
typed_beta <- read.table(typed_beta_path, header=TRUE)
head(typed_beta)     
         SNP Chr     BP A1 A2     BETA     SE NMISS
1 rs12564807   1 734462  G  A  0.95550 0.5353   100
2  rs3094315   1 752566  G  A  0.69770 0.5147   100
3  rs3131972   1 752721  A  G  0.55130 0.5280   100
4  rs3131969   1 754182  A  G -0.04744 0.5306   100
5  rs1048488   1 760912  C  T -0.07650 0.5416   100

# Impute BETA and SE
xx <- SImputeB(ref.geno=ref.geno, ref.map=ref.map, typed=typed_beta, w=1000000, threads=1)

# As discussed above in Zscore imputation, the individual genotype of summary statistics could be used in imputation
xx <- SImputeB(ref.geno=ref.geno, ref.map=ref.map, typed=typed_beta, typed.geno=typed.geno, w=1000000, threads=1)
```

Estimate h2
-----
The estimated heritability (h2) is defined as the coefficient value of regressing square of zscore on ldscore. Compared with REML-based algorithm, LD regression only requires summary statistics and ldscore calculated from reference panel, no need for individual level genoetype data.
```r
sumstat1_path <- system.file("extdata", "sumstat1", package = "SumTool")
ldscore_path <- system.file("extdata", "ldscore", package = "SumTool")
sumstat1 <- read.table(sumstat1_path, header=TRUE)
head(sumstat1)
          SNP Chr   Pos A1 A2       BETA        SE   N
1  rs58108140   1 10583  G  A -0.3328438 0.5004275 100
2 rs189107123   1 10611  C  G -0.3709890 1.2545152 100
3 rs180734498   1 13302  C  T  0.1096603 0.5442463 100
4 rs144762171   1 13327  G  C -0.7596956 0.9316527 100
5 rs151276478   1 13980  T  C -1.0898181 1.2143469 100

ldscore <- read.table(ldscore_path, header=TRUE)
head(ldscore)
          SNP Chr   Pos A1 A2        Maf  ldscore
1  rs58108140   1 10583  G  A 0.20712401 1.203392
2 rs189107123   1 10611  C  G 0.02110818 3.028629
3 rs180734498   1 13302  C  T 0.13720317 1.828657
4 rs144762171   1 13327  G  C 0.03957784 2.650524
5 rs151276478   1 13980  T  C 0.02242744 2.249325

res1 <- LDreg(sumstat = sumstat1, ldscore = ldscore)
```
```
**************************************************
* Summary statistics analysis Tool (SumTool)     *
* Version 0.99.5                                 *
* Author: Lilin Yin                              *
* GPL-3.0 License                                *
**************************************************
Analysis started: 2020-02-25 19:24:37
Number of SNPs for summary statistics 1858
Number of SNPs for ldscore 1851
After merging with reference panel 1851 SNPs remain
Total 664 SNPs that MAF > 0.05
Total 100 individuals included
Using two-step estimator with cutoff at 80
Estimated h2: 0.2052172 (0.03288962)
Lambda GC: 0.8296245
Mean Chi^2: 1.025956
Intercept: 0.4435941 (0.05108769)
Analysis finished: 2020-02-25 19:24:37
Total Running time: 0s
```
***Visualize LD regression results***
```r

```

Estimate rG
-----
In addition to estimate heritability for single trait, the 'LDreg' function can also be applied to estimate the genetic correlation for multiple traits. In this case, just assign a list containing summary statistics of multiple traits to the parameter 'sumstat', the procedure will firstly process the estimation of heritability for each trait and then estimate the genetic correlation for pairs of traits.
```r
sumstat1_path <- system.file("extdata", "sumstat1", package = "SumTool")
sumstat2_path <- system.file("extdata", "sumstat2", package = "SumTool")
ldscore_path <- system.file("extdata", "ldscore", package = "SumTool")
sumstat1 <- read.table(sumstat1_path, header=TRUE)
sumstat2 <- read.table(sumstat2_path, header=TRUE)
ldscore <- read.table(ldscore_path, header=TRUE)

res2 <- LDreg(sumstat = list(sumstat1, sumstat2), ldscore = ldscore)
```
```
**************************************************
* Summary statistics analysis Tool (SumTool)     *
* Version 0.99.5                                 *
* Author: Lilin Yin                              *
* GPL-3.0 License                                *
**************************************************
Analysis started: 2020-02-25 19:26:42
Bivariate LD regression on 2 traits

Heritability of phenotype 1
--------------------------------------------
Number of SNPs for summary statistics 1858
Number of SNPs for ldscore 1851
After merging with reference panel 1851 SNPs remain
Total 664 SNPs that MAF > 0.05
Total 100 individuals included
Using two-step estimator with cutoff at 80
Estimated h2: 0.2052172 (0.03288962)
Lambda GC: 0.8296245
Mean Chi^2: 1.025956
Intercept: 0.4435941 (0.05108769)

Heritability of phenotype 2
--------------------------------------------
Number of SNPs for summary statistics 1858
Number of SNPs for ldscore 1851
After merging with reference panel 1851 SNPs remain
Total 664 SNPs that MAF > 0.05
Total 100 individuals included
Using two-step estimator with cutoff at 80
Estimated h2: 0.2033435 (0.03280185)
Lambda GC: 0.785571
Mean Chi^2: 1.042492
Intercept: 0.451278 (0.05230535)

Genetic Covariance between phenotype 1 and 2
--------------------------------------------
Number of shared SNPs: 1851
Number of individuals for two traits: 100 100
Estimated gencov: 0.1769602 (0.0233377)
Mean z1*z2: 1.024773
Intercept: 0.4565242 (0.04693957)
Genetic Correlation:
          trait1    trait2
trait1 1.0000000 0.8662705
trait2 0.8662705 1.0000000
Analysis finished: 2020-02-25 19:26:42
Total Running time: 0s
```

Estimate Joint Effect
-----
To estimate joint effect, in addition to providing the summary statistics with beta included, the ridge regression lambda also should be calculated in prior. lambda = m * (1 / h2 - 1) where m is the total number of SNPs used in this analysis and h2 (known as heritability) is the proportion of variance in the phenotype explained by all SNPs.
```r
sumstat_path <- system.file("extdata", "typed.marginal", package = "SumTool")
ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")

# load data
sumstat <- read.table(sumstat_path, header=TRUE)
head(sumstat)
         SNP Chr     BP A1 A2     BETA     SE NMISS
1 rs12564807   1 734462  G  A  0.95550 0.5353   100
2  rs3094315   1 752566  G  A  0.69770 0.5147   100
3  rs3131972   1 752721  A  G  0.55130 0.5280   100
4  rs3131969   1 754182  A  G -0.04744 0.5306   100
5  rs1048488   1 760912  C  T -0.07650 0.5416   100
6 rs12562034   1 768448  G  A -0.12360 0.5500   100
data <- read_plink(bfile=ref_bfile_path, threads=1)
geno <- data$geno
map <- data$map
h2 <- 0.5
lambda = nrow(sumstat)*(1/h2-1)
eff <- SBLUP(sumstat = sumstat, geno = geno, map = map, lambda = lambda, w = 1000000, threads = 1)
```
```
**************************************************
* Summary statistics analysis Tool (SumTool)     *
* Version 0.99.5                                 *
* Author: Lilin Yin                              *
* GPL-3.0 License                                *
**************************************************
Analysis started: 2020-02-25 19:30:24
Data and parameters check...(Qualified)
Number of total SNPs in geno is 1851
Number of individuals in geno is 379
Number of total SNPs in summary statistics is 52
Number of individuals in summary statistics is 100
Ridge regression coefficient 52
Number of shared typed SNPs is 45
Window size is 1Mb
Loop on chromosome 1 with 45 SNPs
Total Number of windows: 1
The 1th window: Start[734462] ~ End[998395]
Analysis finished: 2020-02-25 19:30:24
Total Running time: 0s
```
