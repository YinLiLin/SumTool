# SumTool [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/YinLiLin/SumTool/issues) [![](https://img.shields.io/badge/Release-v1.0.0-darkred.svg)](https://github.com/YinLiLin/SumTool) <a href="https://hits.seeyoufarm.com"/><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FYinLiLin%2FSumTool"/></a>

## *A memory-efficient, parallel-accelerated tool for Summary data based Analysis*

Overview
-----
```SumTool``` was designed to implement GWAS summary statistics analysis for big data. It can do ld, ldscore computation, ld pruning and clumping, linear model association, Zscore, Marginal effect imputation, summary data based SNP BLUP (SBLUP) solution, as well as LD regression for heritability and genetic correlation estimation. It is featured with memory efficiency and parallel computation. By the aid of package [bigmemory](https://cran.r-project.org/web/packages/bigmemory/bigmemory.pdf), ```SumTool``` constructs memory-mapped files for genotype panel on disk instead of loading it all into Random Access Memory (RAM), which makes it possible to handle very big data with limited computation resources. Additionally, all matrix manipulations are enhanced by package [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/), which could be sped up by intel MKL library on platform of [Microsoft R Open](https://mran.microsoft.com) (I suggest using MRO rather than R), and all loop procedures are accelerated by OpenMP technology. The functions of SumTool will keep on being enriched with more features based on users feedbacks.

<!--
If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/YinLiLin/SumTool/issues/new):point_left:.
-->
```SumTool``` was developed by [Lilin Yin](https://github.com/YinLiLin) with the support of [Jian Zeng](https://scholar.google.com/citations?user=mOyykToAAAAJ&hl=en) and [Jian Yang](https://scholar.google.com.au/citations?user=aLuqQs8AAAAJ&hl=en). If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/YinLiLin/SumTool/issues/new):point_left:.

Features
-----
- LD 
  - [Linkage disequilibrium matrix](#linkage-disequilibrium)
  - [Linkage disequilibrium score](#ld-score)
  - [Linkage disequilibrium pruning](#ld-pruning)/[clumping](#ld-clumping)
- Association
  - [Linear Model](#linear-model) 
- Imputation
  - [Impute Zscore](#impute-zscore)
  - [Impute Beta and SE](#impute-marginal-effect)
- LD regression
  - [Heritability](#estimate-h2)
  - [Genetic correlation](#estimate-rg)<img src="https://raw.githubusercontent.com/YinLiLin/SumTool/master/pic/P1.png" height="100" align="right" />
- SBLUP
  - [Joint effect](#estimate-joint-effect)
--- 

Installation
-----
Please install ```devtools``` prior to installing ```SumTool```:
```r
> devtools::install_github("YinLiLin/SumTool")
```
After installed successfully, type ```library(SumTool)``` to use, then type ```lsf.str("package:SumTool")``` to see details of all available functions and parameters.

Genotype Converting
-----
With the developing of advaced sequence technologies, the unprecedentedly increased markers across genome make a big challenge in relevant genetic analysis. Loading the genotype into memory directly is highly limited by the computation resources, which becomes the bottleneck of the most of softwares or pipelines. Here we provided two functions ```read_binary, read_vcf``` developed by aid of bigmemory package to construct memory-mapped files ('big.matrix') on disk from plink binary and VCF files. Instead of reading all genotype data into RAM, it greatly reduce memory cost without significantly increase of computation time. Moreover, this data transformation process is only required at the first time, and no matter how big the data is, it can be connected with RAM within seconds at the next time.
```r
# plink binary file
ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
data <- read_binary(bfile=ref_bfile_path, additive=TRUE, threads=1)
# bfile: the prefix of binary files

# or load VCF file
ref_vfile_path <- system.file("extdata", "ref_geno.vcf", package = "SumTool")
data <- read_vcf(vfile=ref_vfile_path, additive=TRUE, threads=1)
# vfile: the full name of vcf file

ref.geno <- data$geno
ref.map <- data$map
```
Missing genotype will be replaced by the major genotype of each allele. By default, the memory-mapped files are directed into work directory, users could redirect to new path as following:
```r
data <- read_binary(bfile=ref_bfile_path, out="./test", threads=1)

# directly use for the next time
ref.geno <- attach.big.matrix("./test.desc")
```
Linkage Disequilibrium
-----
Here, we simply calculate LD by the Pearson correlation (r) of pairs of SNPs
```r
ld <- LDcal(geno=ref.geno, threads=1)
# compute LD for a subset of SNPs
index <- 1:10
geno_sub <- deepcopy(ref.geno, cols=index)
ld_sub <- LDcal(geno_sub, threads=1)
```
By default, it returns a standard square R matrix with dimension of m by m (m is the number of SNPs), which may cost huge space in memory with the increasing SNPs, users can direct a file on disk by parameter 'out' to store the LD matrix in 'big.matrix' format.

LD Score
-----
LD score is defined as the sum of LD r2 between a variant and all the variants in a region. 
```r
ldscore <- LDscore(geno = ref.geno, map = ref.map, w = 100000, b=50000, threads = 1)
```
In ```LDscore```, we provide a parameter 'r2', users could determine to calculate r or r2. By default, the 'r2' is adjusted by r2adj = r2 - [(1 - r2) / (n -2)], as well as 'r'.

LD Pruning
-----
LD pruning is widely used to reduce linked SNPs on base of LD, when it finds a large correlation, it removes one SNP from the correlated pair, keeping the one with the largest minor allele frequency (MAF), thus possibly removing the first SNP. Then it goes on with the next SNP.
```r
snp <- LDprune(geno = ref.geno, map = ref.map, w = 100000, b=50000, threads = 1)
```

LD Clumping
-----
Different with pruning, clumping uses some statistic (usually p-value in the case of GWAS/PRS) to sort the SNPs by importance (e.g. keeping the most significant ones). It takes the first one (e.g. most significant SNP) and removes SNPs that are too correlated with this one in a window around it. As opposed to pruning, this procedure makes sure that this SNP is never removed. We could say that clumping is a trait-specific version of pruning, thus clumping is preferred in Polygenic Risk Score analysis.
```r
p_path <- system.file("extdata", "P.txt", package = "SumTool")
pdata <- read.table(p_path, header = TRUE)
head(pdata)
          SNP         P
1  rs58108140 0.5075209
2 rs189107123 0.7680610
3 rs180734498 0.8407287
4 rs144762171 0.4167849
5 rs151276478 0.3716551
snp <- LDclump(geno = ref.geno, map = ref.map, p = pdata, p.cutoff = 1, r2.cutoff = 0.25, w = 100000, threads = 1)
```

Linear Model
-----
```r
ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
data <- read_binary(bfile=ref_bfile_path, threads=1)
geno <- data$geno
map <- data$map
y <- data$pheno[,1]
gwas <- LMreg(y=y, geno=geno, map=map, threads=1, verbose=FALSE)
```

Impute Zscore
-----
For Zscore imputation, at least 6 columns should be provided in same order with the example for typed SNPs. No need to separate genome into chromosomes. The reference panel can be loaded from both plink binary file (read_binary) or VCF file (read_vcf), we recommend using phased VCF file. Please note that "(..., [additive=FALSE]())" should be added when reading either plink file or VCF file.
```r
# get the path of attached example data on SumTool
ref_file_path <- system.file("extdata", "ref_geno.vcf", package = "SumTool")
typed_z_path <- system.file("extdata", "typed.zscore", package = "SumTool")

# reading data
data <- read_vcf(vfile=ref_file_path, additive=FALSE, threads=1)
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
gwas <- read_binary(bfile=gwas_bfile_path, additive=FALSE, threads=1)
typed.geno <- gwas$geno
typed.map <- gwas$map
# NOTE: the order of SNPs in 'typed.geno' should be consistent with the order in 'typed_z'.
typed.geno <- deepcopy(typed.geno, cols = match(typed_z[, 1], typed.map[, 1]))
xx <- SImputeZ(ref.geno=ref.geno, ref.map=ref.map, typed=typed_z, typed.geno=typed.geno, w=1000000, threads=1)
```

Impute Marginal Effect
-----
For BETA and SE imputation, limited 8 columns should be provided in same order with the example for typed SNPs. No need to separate genome into chromosomes. The reference panel can be loaded from both plink binary file (read_binary) or VCF file (read_vcf), we recommend using phased VCF file. Please note that "(..., [additive=FALSE]())" should be added when reading either plink file or VCF file.<br>
***NOTE***: It is not supported to impute multiple traits at a time.
```r
# get the path of attached example data on SumTool
ref_file_path <- system.file("extdata", "ref_geno.vcf", package = "SumTool")
typed_beta_path <- system.file("extdata", "typed.beta", package = "SumTool")

# reading data
data <- read_vcf(vfile=ref_file_path, additive=FALSE, threads=1)
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
chi2 <- (sumstat1$BETA / sumstat1$SE)^2
ld <- ldscore[match(sumstat1[, 1], ldscore[, 1]), ncol(ldscore)]
ld <- ld[chi2 < 80]
chi2 <- chi2[chi2 < 80]
group <- cut(ld, 50, labels=F)
ld_block <- tapply(ld, group, mean)
chi2_block <- tapply(chi2, group, mean)
plot(ld_block, chi2_block, pch=19, xlab="LD Score Bin", ylab="Mean x^2")
abline(intercept, h2, col="red")
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
To estimate joint effect, in addition to providing the summary statistics with beta included, the ridge regression lambda should be calculated in prior. lambda = m * (1 / h2 - 1) where m is the total number of SNPs used in this analysis and h2 (known as heritability) is the proportion of variance in the phenotype explained by all SNPs.
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
data <- read_binary(bfile=ref_bfile_path, threads=1)
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

Citation
-----
***LD score and regression:***<br>
Bulik-Sullivan, Brendan K., et al. [LD Score regression distinguishes confounding from polygenicity in genome-wide association studies.](https://www.nature.com/articles/ng.3211) ***Nature genetics*** 47.3 (2015): 291.

***Z-score imputation***<br>
Pasaniuc, Bogdan, et al. [Fast and accurate imputation of summary statistics enhances evidence of functional enrichment.](https://academic.oup.com/bioinformatics/article/30/20/2906/2422225) ***Bioinformatics*** 30.20 (2014): 2906-2914.

***Marginal effect imputation***<br>
Zhu, Zhihong, et al. [Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets.](https://www.nature.com/ng/journal/v48/n5/abs/ng.3538.html) ***Nature genetics*** 48.5 (2016): 481.

***SBLUP***<br>
Robinson, Matthew R., et al. [Genetic evidence of assortative mating in humans.](https://www.nature.com/articles/s41562-016-0016?fbclid=IwAR1hYuVgzjT1S4dMixO9-_IHRgCoFncZBBEU5ASCly0cN_9Kmk4ZtLGY3fY) ***Nature Human Behaviour*** 1.1 (2017): 0016.

