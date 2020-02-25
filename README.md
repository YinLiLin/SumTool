# SumTool [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/YinLiLin/SumTool/issues) [![](https://img.shields.io/badge/Release-v0.99.5-darkred.svg)](https://github.com/YinLiLin/SumTool)

## *A memory-efficient, parallel-accelerated tool for Summary data based Analysis*

Overview
-----
```SumTool``` was designed to do analysis using GWAS summary statistics for big data. It can do ld, ldscore computation, Zscore, Marginal effect imputation, summary data based SNP BLUP (SBLUP) solution, as well as LD regression for heritability and genetic correlation estimation. It featured with memory efficiency and parallel calculation. By the aid of package 'bigmemory', SumTool constructs memory-mapped files for genotype panel on disk instead of loading it all into Random Access Memory (RAM), which makes it possible to handle very big data with limited computation resources. Additionally, all parallel processes are accelerated as fast as possible by OpenMP technology. The functions of SumTool will keep on being enriched with more features based on users feedbacks.

```SumTool``` was developed by [Lilin Yin](https://github.com/YinLiLin) with the support of [Jian Zeng](https://scholar.google.com/citations?user=mOyykToAAAAJ&hl=en) and [Jian Yang](https://scholar.google.com.au/citations?user=aLuqQs8AAAAJ&hl=en). If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/YinLiLin/SumTool/issues/new):point_left:.

Contents
-----

- LD 
  - [Linkage disequilibrium matrix](#linkage-disequilibrium)
  - [Linkage disequilibrium score](#ld-score)
- Impute
  - [Zscore imputation](#impute-zscore)
  - [Beta and SE imputation](#impute-marginal-effect)
- LD regression
  - [Heritability](#estimate-h2)
  - [Genetic correlation](#estimate-rg)
- SBLUP
  - [Joint effect](#estimate-joint-effect)
  
  
  
Linkage Disequilibrium
-----




LD Score
-----




Impute Zscore
-----



Impute Marginal Effect
-----


Impute Zscore
-----



Estimate h2
-----



Estimate rG
-----




Estimate Joint Effect
-----


