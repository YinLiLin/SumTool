version.info <- function(width=50)
{
	if (getPackageName() == ".GlobalEnv") {
        version <- "devel"
    } else {
        version <- as.character(packageVersion(getPackageName()))
    }
    Author <- "Lilin Yin [ylilin@163.com]"
	cat(rep("*", width), "\n", sep="")
	cat("* Summary statistic analysis Tool (SumTool)", rep(" ", width-44), "*", "\n", sep="")
	cat("* Version ", version, rep(" ", width-nchar(version)-11), "*", "\n", sep="")
	cat("* Author: ", Author, rep(" ", width-11-nchar(Author)), "*", "\n", sep="")
	cat("* GPL-3.0 License", rep(" ", width-18), "*", "\n", sep="")
	cat(rep("*", width), "\n", sep="")
}

times <-
function(x)
{
	h <- x %/% 3600
	m <- (x %% 3600) %/% 60
	s <- ((x %% 3600) %% 60)
	index <- which(c(h, m, s) != 0)
	num <- c(h, m, s)[index]
	num <- round(num, 0)
	char <- c("h", "m", "s")[index]
	return(paste(num, char, sep="", collapse=""))
}

#' Data reading
#'
#' To read Plink binary format data into bigmemory (0, 1, 2) format
#'
#' @param bfile character, prefix of Plink binary format data.
#' @param maxLine number, set the number of lines to handle at a time, bigger lines require more memory.
#' @param impute logical, whether to impute missing values in genotype.
#' @param backingpath the path to the directory containing the file backing cache.
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with ‘attach.big.matrix’; if ‘NULL’, the ‘backingfile’ is used as the root part of the descriptor file name.  The descriptor file is placed in the same directory as the backing files.
#' @param backingfile the root name for the file(s) for the cache.
#' @param verbose logical, whether to print the log information
#' @param threads number, the number of used threads for parallel process

#' @examples
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' geno <- data$geno
#' map <- data$map

read_plink <- 
function(
	bfile = "", 
	maxLine = 1000,
	impute = TRUE,
	backingpath = NULL,
	descriptorfile = NULL,
	backingfile = NULL,
	verbose = TRUE,
	threads = 1
){
	if(verbose)	version.info()
	t1 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	m <- FileNrow(paste(bfile, ".bim", sep=""))
	n <- FileNrow(paste(bfile, ".fam", sep=""))
	if(verbose)	cat("Number of SNPs: ", m , "\n", sep="")
	if(verbose)	cat("Number of individuals: ", n , "\n", sep="")
	if(verbose)	cat("Reading...\n")
	if(!is.null(backingfile))	backingfile <- basename(backingfile)
	if(!is.null(descriptorfile))	descriptorfile <- basename(descriptorfile)
	if(!is.null(backingfile) & is.null(descriptorfile))
		stop("descriptorfile shoud be assigned.")
	if(is.null(backingfile) & !is.null(descriptorfile))
		stop("backingfile shoud be assigned.")
	map_file = NULL
	if(!is.null(descriptorfile)){
		map_file <- unlist(strsplit(descriptorfile, "", fixed = TRUE))
		sep_index <- which(map_file == ".")
		if(length(sep_index)){
			map_file <- paste0(map_file[1 : (sep_index[length(sep_index)] - 1)], collapse="")
		}else{
			map_file <- paste0(map_file, collapse="")
		}
		if(!is.null(backingpath))	map_file <- paste0(backingpath, "/", map_file)
	}
	map <- as.data.frame(rMap_c(paste0(bfile, ".bim"), out = map_file), stringsAsFactors=FALSE)
	map$Pos <- as.numeric(map$Pos)
	pheno <- read.table(paste0(bfile, ".fam"), header=FALSE)[, -c(1:5), drop=FALSE]
	geno <- bigmemory::big.matrix(
		nrow = n,
		ncol = m,
		type = "char",
		dimnames = c(NULL, NULL),
		backingfile = backingfile,
		backingpath = backingpath,
		descriptorfile = descriptorfile
	)
	rData_c(bfile = bfile, pBigMat = geno@address, maxLine = maxLine, impt = impute, verbose = verbose, threads = threads)
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2 - t1), "\n")
	return(list(pheno=pheno, geno=geno, map=map))	
}

#' Data writing
#'
#' To write bigmemory (0, 1, 2) data into Plink binary format
#'
#' @param geno big.matrix, numeric matrix stored in big.matrix format.
#' @param map number, set the number of lines to handle at a time, bigger lines require more memory.
#' @param out logical, whether to impute missing values in genotype.
#' @param verbose logical, whether to print the log information
#' @param threads number, the number of used threads for parallel process

#' @examples
#' # reading data
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' geno <- data$geno
#' map <- data$map
#'
#' # writing data
#' # write_plink(geno=geno, map=map, out="./test", threads=1)

write_plink <- 
function(
	geno = NULL, 
	map = NULL,
	out = "plink",
	verbose = TRUE,
	threads = 1
){
	if(verbose)	version.info()
	t1 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	if(is.null(geno))	stop("genotype in big.matrix should be provided!")
	if(!is.big.matrix(geno))	stop("genotype should be in 'big.matrix' format!")
	if(is.null(map))	stop("map information should be provided!")
	if(ncol(map) != 5)	stop("5 columns should be provided in map!")
	if(ncol(geno) != nrow(map))	stop("number of SNPs not equals between genotype and map!")
	m <- ncol(geno)
	n <- nrow(geno)
	if(verbose)	cat("Number of SNPs: ", m , "\n", sep="")
	if(verbose)	cat("Number of individuals: ", n , "\n", sep="")
	if(verbose)	cat("Writing...\n")
	wData_c(pBigMat = geno@address, bed_file = out, verbose = verbose, threads = threads)
	write.table(cbind(1 : n, 1 : n, 0, 0, 0, -9), paste0(out, '.fam'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ' ')
	write.table(cbind(map[, 2], map[, 1], 0, map[, 3 : 5]), paste0(out, '.bim'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ' ')
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2 - t1), "\n")
	return(NULL)	
}

#' Summary data based imputation on Z-score
#'
#' To impute Zscore using summary data by SImpute/SImpute-LD
#'
#' @param ref.geno big.matrix (n1 * m1), reference genotype panel
#' @param ref.map matrix (m1 * 5): SNPs, Chr, position, A1, A2
#' @param typed.geno big.matrix (n2 * m2), individual level genotype for typed SNPs (this file is optional)
#' @param typed matrix (m2 * 6): SNPs, Chr, position, A1, A2, Z
#' @param w number, set the window size in bp, default 1000000
#' @param b number, set the buffer size in bp of each window, default 250000
#' @param lambda number, ridge regression value on LD matrix of typed SNPs: solve(Rtt + diag(lambda))
#' @param maf number, SNPs whose minor allele frequency are lower than set value will not be imputed
#' @param correlation logical, if TRUE, the LD matrix will be constructed by correlation of all pairs
#' @param verbose logical, whether to print the log information
#' @param threads number, the number of used threads for parallel process

#' @examples
#' #----------------SImpute----------------#
#'
#' # get path of the attached files in package
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' typed_z_path <- system.file("extdata", "typed.zscore", package = "SumTool")
#'
#' # reading data
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' ref.geno <- data$geno
#' ref.map <- data$map
#' typed_z <- read.table(typed_z_path, header=TRUE)
#'
#' # Impute Zscore
#' xx <- SImputeZ(ref.geno=ref.geno, ref.map=ref.map, typed=typed_z, threads=1)
#'
#' #--------------SImpute-LD---------------#
#'
#' gwas_bfile_path <- system.file("extdata", "gwas_geno", package = "SumTool")
#' gwas <- read_plink(bfile=gwas_bfile_path, threads=1)
#' typed.geno <- gwas$geno
#' typed.map <- gwas$map
#' #NOTE: the order of SNPs in 'typed.geno' should be consistent with the order in 'typed_z'.
#' typed.geno <- deepcopy(typed.geno, cols = match(typed_z[, 1], typed.map[, 1]))
#'
#' # Impute Zscore
#' xx <- SImputeZ(ref.geno=ref.geno, ref.map=ref.map, typed=typed_z, typed.geno=typed.geno, threads=1)

SImputeZ <- function(ref.geno = NULL, ref.map = NULL, typed.geno = NULL, typed = NULL, w = 1000000, b = 500000, lambda = NULL, maf = 0.000001, correlation = TRUE, verbose = TRUE, threads = 1)
{

	SImpute_bin <- function(ref.geno, ref.map, typed.geno = NULL, typed, lambda = 0.001, maf = 0.000001, correlation = TRUE, verbose = TRUE, threads = 1)
	{
		index <- match(typed[, 1], ref.map[, 1])
		# if(sum(is.na(index)) != 0)	stop("Some typed SNPs don't exit in reference genotype panel!")
		# typed <- typed[order(index), ]
		# if(!is.null(typed.geno))	typed.geno <- deepcopy(typed.geno, cols = order(index))
		# index <- match(typed[, 1], ref.map[, 1])
		if(is.null(typed.geno)){
			imp <- SImputeZ_bin_c(ref.geno@address, typed_index = index, typed_value = data.matrix(typed[, -c(1:5), drop=FALSE]), verbose = verbose, lambda = lambda, maf = maf, haps = !correlation, threads = threads)
		}else{
			imp <- SImputeZ_ld_bin_c(ref.geno@address, typed.geno@address, typed_index = index, typed_value = data.matrix(typed[, -c(1:5), drop=FALSE]), verbose = verbose, lambda = lambda, maf = maf, haps = !correlation, threads = threads)
		}
		imp <- data.frame(ref.map, imp)
	}

	if(verbose)	version.info()
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	#check parameter
	if(verbose)	cat("Data and parameters check...")
	if(is.null(ref.geno))	stop("Please provide ref.geno!")
	if(!is.big.matrix(ref.geno))	stop("genotype should be in 'big.matrix' format!")
	if(is.null(ref.map))	stop("Please provide ref.map!")
	if(is.null(typed))	stop("Please provide typed!")
	if(any(duplicated(ref.map[, 1])))	stop("duplicated SNP names exist in ref.map file!")
	if(any(duplicated(typed[, 1])))	stop("duplicated SNP names exist in typed file!")
	if(ncol(ref.map) != 5)	stop("Only 5 columns limited for ref.map! (snp, chr, pos, a1, a2)")
	if(ncol(typed) < 6)	stop("At least 6 columns should be provided for typed! (snp, chr, pos, a1, a2, zscore")
	if(!is.numeric(ref.map[, 3]))	stop("Physical position in map should be numeric defined!")
	if(!is.numeric(typed[, 3]))	stop("Physical position in typed should be numeric defined!")
	for(i in 1 : ncol(ref.map)){
		if(is.factor(ref.map[, i]))	ref.map[, i] <- as.character.factor(ref.map[, i])
		if(is.factor(typed[, i]))	typed[, i] <- as.character.factor(typed[, i])
	}	
	if(sum(is.na(typed[, -c(1:5)])) != 0)	stop(paste("NA is not allowed in Zscore!", sep=""))
	if(ncol(ref.geno) != nrow(ref.map))	stop("Number of SNPs not equals between ref.geno and ref.map!")
	#if(hasNA(ref.geno@address, threads = threads))	stop("NA is not allowed in ref.geno!")
	if(!is.null(typed.geno)){
		if(!is.big.matrix(typed.geno))	stop("genotype should be in 'big.matrix' format!")
		if(ncol(typed.geno) != nrow(typed))	stop("Number of SNPs not equals between typed.geno and typed!")
	}
	#if(!is.null(typed.geno) && hasNA(typed.geno@address, threads = threads))	stop("NA is not allowed in typed.geno!")

	SNP_NA <- is.na(ref.map[, 3]) | ref.map[, 3] == 0 | ref.map[, 3] == -9
	if(sum(SNP_NA) != 0){
		ref.geno <- deepcopy(ref.geno, cols = !SNP_NA)
		ref.map <- ref.map[!SNP_NA, ]
	}
	typed_N <- nrow(typed)
	SNP_NA <- is.na(typed[, 3]) | typed[, 3] == 0 | typed[, 3] == -9
	if(sum(SNP_NA) != 0){
		if(!is.null(typed.geno)) typed.geno <- deepcopy(typed.geno, cols = !SNP_NA)
		typed <- typed[!SNP_NA, , drop=FALSE]
	}
	mergedSNP <- intersect(typed[, 1], ref.map[, 1])
	if(length(mergedSNP) == 0)	stop("No shared SNPs between reference and typed SNPs!")
	if(length(mergedSNP) != nrow(typed)){
		indx <- typed[, 1] %in% mergedSNP
		typed_temp <- typed[!indx, ]
		colnames(typed_temp) <- colnames(typed)
		typed <- typed[indx, ]
		if(!is.null(typed.geno)) typed.geno <- deepcopy(typed.geno, cols = indx)
	}else{
		typed_temp <- NULL
	}
	row_logi <- NULL
	for(i in 2:ncol(ref.map)){
		row_logi <- cbind(row_logi, ref.map[match(mergedSNP, ref.map[, 1]), i] == typed[, i])
	}
	typed_index <- apply(row_logi, 1, all); rm(row_logi); gc()
	# typed_index <- sapply(1:nrow(typed), function(i){
	# 	if(all(typed[i, 1:4] == ref.map[typed_match_index[i], ])){
	# 		return(TRUE)
	# 	}else{
	# 		return(FALSE)
	# 	}
	# })

	if(verbose) cat("(Qualified)\n")
	if(!all(typed_index)){
		typed <- typed[typed_index, , drop=FALSE]
		if(!is.null(typed.geno)) typed.geno <- deepcopy(typed.geno, cols = typed_index)
		if(verbose)	cat("(Warning:", sum(!typed_index), "SNPs sharing same names between reference and typed don't have equal map information!)\n")
	}

	#SImpute or SImpute-LD
	if(!is.null(typed.geno)){
		if(is.null(lambda))	lambda <- 0.001
		if(verbose)	cat(paste("SImpute-LD on Zscore started...\n", sep=""))
	}else{
		if(is.null(lambda))	lambda <- 0.1
		if(verbose)	cat(paste("SImpute on Zscore started...\n", sep=""))
	}
	gc()

	t1 <- as.numeric(Sys.time())

	#confirm parameters
	if(verbose)	cat("Number of traits is ", ncol(typed) - 5, "\n", sep="")
	if(verbose)	cat("Number of total SNPs in reference is ", ncol(ref.geno) , "\n", sep="")
	if(verbose)	cat("Number of individuals in reference is ", nrow(ref.geno) , "\n", sep="")
	if(verbose)	cat("Number of typed SNPs is ", typed_N, "\n", sep="")
	if(verbose)	cat("Number of shared typed SNPs is ", nrow(typed), "\n", sep="")
	if(verbose && !is.null(typed.geno))	cat("Number of individuals in GWAS sample ", nrow(typed.geno), "\n", sep="")
	if(verbose)	cat("MAF threshold is ", maf, "\n", sep="")
	if(verbose)	cat("Window size is ", w / 1e6, "Mb", "\n", sep="")
	if(verbose)	cat("Buffer size is ", b / 1e3, "Kb", "\n", sep="")

	if(verbose){
		if(!correlation)	cat("Using haplotype for LD matrix\n")
		if(correlation)	cat("Using correlation for LD matrix\n")
	}
	if(verbose)	cat("Lambda in ridge regression is", lambda, "\n")

	chr <- unique(ref.map[, 2])

	imp <- typed_temp
	for(chri in chr){
		chri_index <- ref.map[, 2] == chri
		if(verbose)	cat("Loop on chromosome", chri, "with", sum(chri_index), "SNPs\n")
		typed_chri_index <- typed[, 2] == chri
		# chri_pos_min <- min(ref.map[chri_index, 3])
		chri_pos_min <- 1
		chri_pos_max <- max(ref.map[chri_index, 3])
		loop <- TRUE
		wind_min <- 1
		wind_max <- w
		# if(verbose)	cat("Divide genome into windows...\n")
		while(loop){
			wind_min <- c(wind_min, wind_max[length(wind_max)] + 1)
			wind_max <- c(wind_max, wind_max[length(wind_max)] + w)
			if((wind_max[length(wind_max)] + b / 2) >= chri_pos_max)	loop <- FALSE
		}
		if(w >= chri_pos_max){
			wind_min <- wind_min[1]
			wind_max <- wind_max[1]
		}
		if(verbose)	cat("Total Number of windows: ", length(wind_min), "\n", sep="")
		# if(verbose)	cat("Impute for windows...\n")
		wind_max[length(wind_max)] <- chri_pos_max
		wind_min_b <- wind_min - b / 2
		wind_max_b <- wind_max + b / 2
		wind_min_b[1] <- 1
		wind_max_b[length(wind_max_b)] <- chri_pos_max
		simpute_mc <- function(i){
			index1 <- which(chri_index & ref.map[, 3] >= wind_min_b[i] & ref.map[, 3] <= wind_max_b[i])
			index2 <- typed_chri_index & typed[, 3] >= wind_min_b[i] & typed[, 3] <= wind_max_b[i]
			if(length(index1) != 0 & sum(index2) != 0){
				if(verbose)	cat("The ", i, "th window of Chr ", chri, ": Start[", min(ref.map[index1, 3], na.rm = TRUE),"] ~ End[", max(ref.map[index1, 3], na.rm = TRUE),"]\n", sep="")
				if(length(index1) == nrow(ref.map)){
					if(is.null(typed.geno)){
						imp_bin <- SImpute_bin(ref.geno = ref.geno, ref.map = ref.map, typed.geno = typed.geno, typed = (typed[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}else{
						imp_bin <- SImpute_bin(ref.geno = ref.geno, ref.map = ref.map, typed.geno = deepcopy(typed.geno, cols = index2), typed = (typed[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}
				}else{
					if(is.null(typed.geno)){
						imp_bin <- SImpute_bin(ref.geno = (deepcopy(ref.geno, cols=index1)), ref.map = (ref.map[index1, ]), typed.geno = typed.geno, typed = (typed[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}else{
						imp_bin <- SImpute_bin(ref.geno = (deepcopy(ref.geno, cols=index1)), ref.map = (ref.map[index1, ]), typed.geno = deepcopy(typed.geno, cols = index2), typed = (typed[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}
				}
			}else if(length(index1) != 0 & sum(index2) == 0){
				imp_bin <- cbind(ref.map[index1, ], matrix(NA, length(index1), (ncol(typed) - 5)))
			}else{
				imp_bin <- NULL
			}
			if(length(wind_min) != 1 && !is.null(imp_bin)){
				index <- (imp_bin[, 3] >= wind_min[i]) & (imp_bin[, 3] <= wind_max[i])
				imp_bin <- imp_bin[index, ]
			}
			if(!is.null(imp_bin))	colnames(imp_bin) <- colnames(typed)
			return(imp_bin)
		}
		imp_bin <- lapply(1 : length(wind_min_b), simpute_mc)
		if(verbose)	cat("Merge imputed value within windows...\n")
	
		imp <- rbind(imp, do.call(rbind, imp_bin))
	}
	imp <- imp[order(imp[, 2], imp[, 3]), ]
	if(verbose)	cat("Imputation accomplished successfully.\n")
	colnames(imp)[1 : 5] <- c("SNP", "Chr", "Pos", "A1", "A2")
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(imp)
}

#' Summary data based imputation on marginal effect
#'
#' To impute marginal effect using summary data by SImpute/SImpute-LD
#'
#' @param ref.geno big.matrix (n1 * m1), reference genotype panel
#' @param ref.map matrix (m1 * 5): SNPs, Chr, position, A1, A2
#' @param typed.geno big.matrix (n2 * m2), individual level genotype for typed SNPs (this file is optional)
#' @param typed matrix (m2 * 8): SNPs, Chr, position, A1, A2, BETA, SE, N
#' @param w number, set the window size in bp, default 1000000
#' @param b number, set the buffer size in bp of each window, default 250000
#' @param lambda number, ridge regression value on LD matrix of typed SNPs: solve(Rtt + diag(lambda))
#' @param maf number, SNPs whose minor allele frequency are lower than set value will not be imputed
#' @param correlation logical, if TRUE, the LD matrix will be constructed by correlation of all pairs
#' @param verbose logical, whether to print the log information
#' @param threads number, the number of used threads for parallel process

#' @examples
#' #----------------SImpute----------------#
#'
#' # get path of the attached files in package
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' typed_b_path <- system.file("extdata", "typed.marginal", package = "SumTool")
#'
#' # reading data
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' ref.geno <- data$geno
#' ref.map <- data$map
#' typed_b <- read.table(typed_b_path, header=TRUE)
#'
#' # Impute marginal effect and se
#' xx <- SImputeB(ref.geno=ref.geno, ref.map=ref.map, typed=typed_b, threads=1)
#'
#' #--------------SImpute-LD---------------#
#'
#' gwas_bfile_path <- system.file("extdata", "gwas_geno", package = "SumTool")
#' gwas <- read_plink(bfile=gwas_bfile_path, threads=1)
#' typed.geno <- gwas$geno
#' typed.map <- gwas$map
#' #NOTE: the order of SNPs in 'typed.geno' should be consistent with the order in 'typed_b'.
#' typed.geno <- deepcopy(typed.geno, cols = match(typed_b[, 1], typed.map[, 1]))
#'
#' # Impute marginal effect and se
#' xx <- SImputeB(ref.geno=ref.geno, ref.map=ref.map, typed=typed_b, typed.geno=typed.geno, threads=1)

SImputeB <- function(ref.geno = NULL, ref.map = NULL, typed.geno = NULL, typed = NULL, w = 1000000, b = 500000, lambda = NULL, maf = 0.000001, correlation = TRUE, verbose = TRUE, threads = 1)
{

	SImpute_bin <- function(ref.geno, ref.map, typed.geno = NULL, typed, lambda = 0.001, maf = 0.000001, correlation = TRUE, verbose = TRUE, threads = 0)
	{
		index <- match(typed[, 1], ref.map[, 1])
		# if(sum(is.na(index)) != 0)	stop("Some typed SNPs don't exit in reference genotype panel!")
		# typed <- typed[order(index), ]
		# if(!is.null(typed.geno))	typed.geno <- deepcopy(typed.geno, cols = order(index))
		# index <- match(typed[, 1], ref.map[, 1])
		if(is.null(typed.geno)){
			imp <- SImputeZ_bin_c(ref.geno@address, typed_index = index, typed_value = data.matrix(typed[, -c(1:5), drop=FALSE]), verbose = verbose, lambda = lambda, maf = maf, haps = !correlation, threads = threads)
		}else{
			imp <- SImputeZ_ld_bin_c(ref.geno@address, typed.geno@address, typed_index = index, typed_value = data.matrix(typed[, -c(1:5), drop=FALSE]), verbose = verbose, lambda = lambda, maf = maf, haps = !correlation, threads = threads)
		}
		imp <- data.frame(ref.map, imp)
	}

	if(verbose)	version.info()
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	
	#check parameter
	if(verbose)	cat("Data and parameters check...")
	if(is.null(ref.geno))	stop("Please provide ref.geno!")
	if(!is.big.matrix(ref.geno))	stop("genotype should be in 'big.matrix' format!")
	if(is.null(ref.map))	stop("Please provide ref.map!")
	if(is.null(typed))	stop("Please provide typed!")
	if(any(duplicated(ref.map[, 1])))	stop("duplicated SNP names exist in ref.map file!")
	if(any(duplicated(typed[, 1])))	stop("duplicated SNP names exist in typed file!")
	if(ncol(ref.map) != 5)	stop("Only 5 columns limited for ref.map! (snp, chr, pos, a1, a2)")
	if(ncol(typed) != 8)	stop("Only 8 columns should be provided for typed! (snp, chr, pos, a1, a2, beta, se, n)")
	if(!is.numeric(ref.map[, 3]))	stop("Physical position in map should be numeric defined!")
	if(!is.numeric(typed[, 3]))	stop("Physical position in typed should be numeric defined!")
	for(i in 1 : ncol(ref.map)){
		if(is.factor(ref.map[, i]))	ref.map[, i] <- as.character.factor(ref.map[, i])
		if(is.factor(typed[, i]))	typed[, i] <- as.character.factor(typed[, i])
	}
	if(sum(is.na(typed[, -c(1:5)])) != 0)	stop(paste("NA is not allowed in typed information!", sep=""))
	if(ncol(ref.geno) != nrow(ref.map))	stop("Number of SNPs not equals between ref.geno and ref.map!")
	#if(hasNA(ref.geno@address, threads = threads))	stop("NA is not allowed in ref.geno!")
	if(!is.null(typed.geno)){
		if(!is.big.matrix(typed.geno))	stop("genotype should be in 'big.matrix' format!")
		if(ncol(typed.geno) != nrow(typed))	stop("Number of SNPs not equals between typed.geno and typed!")
	}
	#if(!is.null(typed.geno) && hasNA(typed.geno@address, threads = threads))	stop("NA is not allowed in typed.geno!")

	SNP_NA <- is.na(ref.map[, 3]) | ref.map[, 3] == 0 | ref.map[, 3] == -9
	if(sum(SNP_NA) != 0){
		ref.geno <- deepcopy(ref.geno, cols = !SNP_NA)
		ref.map <- ref.map[!SNP_NA, ]
	}
	typed_N <- nrow(typed)
	SNP_NA <- is.na(typed[, 3]) | typed[, 3] == 0 | typed[, 3] == -9
	if(sum(SNP_NA) != 0){
		if(!is.null(typed.geno)) typed.geno <- deepcopy(typed.geno, cols = !SNP_NA)
		typed <- typed[!SNP_NA, , drop=FALSE]
	}
	mergedSNP <- intersect(typed[, 1], ref.map[, 1])
	if(length(mergedSNP) == 0)	stop("No shared SNPs between reference and typed SNPs!")
	if(length(mergedSNP) != nrow(typed)){
		indx <- typed[, 1] %in% mergedSNP
		typed_temp <- typed[!indx, -ncol(typed)]
		colnames(typed_temp)[1:5] <- c("SNP", "Chr", "Pos", "A1", "A2")
		typed <- typed[indx, ]
		if(!is.null(typed.geno)) typed.geno <- deepcopy(typed.geno, cols = indx)
	}else{
		typed_temp <- NULL
	}
	row_logi <- NULL
	for(i in 2:ncol(ref.map)){
		row_logi <- cbind(row_logi, ref.map[match(mergedSNP, ref.map[, 1]), i] == typed[, i])
	}
	typed_index <- apply(row_logi, 1, all); rm(row_logi); gc()
	# typed_index <- sapply(1:nrow(typed), function(i){
	# 	if(all(typed[i, 1:4] == ref.map[typed_match_index[i], ])){
	# 		return(TRUE)
	# 	}else{
	# 		return(FALSE)
	# 	}
	# })
	if(verbose) cat("(Qualified)\n")
	if(!all(typed_index)){
		typed <- typed[typed_index, , drop=FALSE]
		if(!is.null(typed.geno)) typed.geno <- deepcopy(typed.geno, cols = typed_index)
		if(verbose)	cat("(Warning:", sum(!typed_index), "SNPs sharing same names between reference and typed don't have equal map information!)\n")
	}

	#SImpute or SImpute-LD
	if(!is.null(typed.geno)){
		n <- nrow(typed.geno)
		if(is.null(lambda))	lambda <- 0.001
		if(verbose)	cat(paste("SImpute-LD on Marginal effect and SE started...\n", sep=""))
	}else{
		n <- round(mean(typed[, ncol(typed)]))
		if(is.null(lambda))	lambda <- 0.1
		if(verbose)	cat(paste("SImpute on Marginal effect and SE started...\n", sep=""))
	}
	gc()
	
	t1 <- as.numeric(Sys.time())

	xixj_all <- BigStat(ref.geno@address, index_ = 0 : (ncol(ref.geno) - 1), threads = threads)$sd
	xixj_all <- xixj_all * xixj_all
	if(!is.null(typed.geno)){
		xixj_typed <- BigStat(typed.geno@address, index_ = 0 : (ncol(typed.geno) - 1), threads = threads)$sd
		xixj_typed <- xixj_typed * xixj_typed
	}else{
		xixj_typed <- xixj_all[match(typed[, 1], ref.map[, 1])] * (n / nrow(ref.geno))
	}

	vary <- mean(xixj_typed * (n * (typed[, 7])^2 + (typed[, 6])^2) / (n - 1))
	z <- typed[, 6] / typed[, 7]

	#confirm parameters
	if(verbose)	cat("Phenotypic variance is", vary, "\n")
	if(verbose)	cat("Number of total SNPs in reference is ", ncol(ref.geno) , "\n", sep="")
	if(verbose)	cat("Number of individuals in reference is ", nrow(ref.geno) , "\n", sep="")
	if(verbose)	cat("Number of typed SNPs is ", typed_N, "\n", sep="")
	if(verbose)	cat("Number of shared typed SNPs is ", nrow(typed), "\n", sep="")
	if(verbose) cat("Number of individuals in summary data ", n, "\n", sep="")
	if(verbose)	cat("MAF threshold is ", maf, "\n", sep="")
	if(verbose)	cat("Window size is ", w / 1e6, "Mb", "\n", sep="")
	if(verbose)	cat("Buffer size is ", b / 1e3, "Kb", "\n", sep="")

	if(verbose){
		if(!correlation)	cat("Using haplotype for LD matrix\n")
		if(correlation)	cat("Using correlation for LD matrix\n")
	}
	if(verbose)	cat("Lambda in ridge regression is", lambda, "\n")
	chr <- unique(ref.map[, 2])
	imp <- typed_temp
	for(chri in chr){
		chri_index <- ref.map[, 2] == chri
		if(verbose)	cat("Loop on chromosome", chri, "with", sum(chri_index), "SNPs\n")
		typed_chri_index <- typed[, 2] == chri
		# chri_pos_min <- min(ref.map[chri_index, 3])
		chri_pos_min <- 1
		chri_pos_max <- max(ref.map[chri_index, 3])
		loop <- TRUE
		wind_min <- 1
		wind_max <- w
		# if(verbose)	cat("Divide genome into windows...\n")
		while(loop){
			wind_min <- c(wind_min, wind_max[length(wind_max)] + 1)
			wind_max <- c(wind_max, wind_max[length(wind_max)] + w)
			if((wind_max[length(wind_max)] + b / 2) >= chri_pos_max)	loop <- FALSE
		}
		if(w >= chri_pos_max){
			wind_min <- wind_min[1]
			wind_max <- wind_max[1]
		}
		if(verbose)	cat("Total Number of windows: ", length(wind_min), "\n", sep="")
		# if(verbose)	cat("Impute for windows...\n")
		wind_max[length(wind_max)] <- chri_pos_max
		wind_min_b <- wind_min - b / 2
		wind_max_b <- wind_max + b / 2
		wind_min_b[1] <- 1
		wind_max_b[length(wind_max_b)] <- chri_pos_max
		simpute_mc <- function(i){
			index1 <- which(chri_index & ref.map[, 3] >= wind_min_b[i] & ref.map[, 3] <= wind_max_b[i])
			index2 <- typed_chri_index & typed[, 3] >= wind_min_b[i] & typed[, 3] <= wind_max_b[i]
			if(length(index1) != 0 & sum(index2) != 0){
				if(verbose)	cat("The ", i, "th window of Chr ", chri, ": Start[", min(ref.map[index1, 3], na.rm = TRUE),"] ~ End[", max(ref.map[index1, 3], na.rm = TRUE),"]\n", sep="")
				if(length(index1) == nrow(ref.map)){
					if(is.null(typed.geno)){
						imp_bin <- SImpute_bin(ref.geno = ref.geno, ref.map = ref.map, typed.geno = typed.geno, typed = (cbind(typed[, 1 : 5], z)[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}else{
						imp_bin <- SImpute_bin(ref.geno = ref.geno, ref.map = ref.map, typed.geno = deepcopy(typed.geno, cols = index2), typed = (cbind(typed[, 1 : 5], z)[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}
				}else{
					if(is.null(typed.geno)){
						imp_bin <- SImpute_bin(ref.geno = (deepcopy(ref.geno, cols=index1)), ref.map = (ref.map[index1, ]), typed.geno = typed.geno, typed = (cbind(typed[, 1 : 5], z)[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}else{
						imp_bin <- SImpute_bin(ref.geno = (deepcopy(ref.geno, cols=index1)), ref.map = (ref.map[index1, ]), typed.geno = deepcopy(typed.geno, cols = index2), typed = (cbind(typed[, 1 : 5], z)[index2, ]), verbose = verbose, lambda = lambda, maf = maf, correlation = correlation, threads = threads)
					}
				}
			}else if(length(index1) != 0 & sum(index2) == 0){
				imp_bin <- cbind(ref.map[index1, ], matrix(NA, length(index1), 1))
			}else{
				imp_bin <- NULL
			}
			if(length(wind_min) != 1 && !is.null(imp_bin)){
				index <- (imp_bin[, 3] >= wind_min[i]) & (imp_bin[, 3] <= wind_max[i])
				imp_bin <- imp_bin[index, ]
			}
			if(!is.null(imp_bin))	colnames(imp_bin) <- c("SNP", "Chr", "Pos", "A1", "A2", "Z")
			return(imp_bin)
		}
		imp_bin <- lapply(1 : length(wind_min_b), simpute_mc)
		if(verbose)	cat("Merge imputed value within windows...\n")
		impz <- do.call(rbind, imp_bin)
		impb <- cbind(impz[, 1 : 5], matrix(NA, nrow(impz), 2))
		scaled <- sqrt((xixj_all[match(impz[, 1], ref.map[, 1])] * n / nrow(ref.geno)) * (n - 1 + (impz[, 6])^2))
		impb[, 6] <- impz[, 6] * sqrt(vary * (n - 1)) / scaled
		impb[, 7] <- sqrt(vary * (n - 1)) / scaled
		colnames(impb) <- c("SNP", "Chr", "Pos", "A1", "A2", colnames(typed)[6:7])
		rm(imp_bin); rm(impz); gc();
		imp <- rbind(imp, impb)
	}
	imp[match(typed[, 1], imp[, 1]), c(6:7)] <- typed[, c(6:7)]
	imp <- imp[order(imp[, 2], imp[, 3]), ]
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Imputation accomplished successfully.\n")
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(imp)
}

#' LD calculation
#'
#' To calculate LD for all pairs of SNPs on given genome
#'
#' @param geno big.matrix (n * m), genotype
#' @param index vector, only calculate LD for a subset with all SNPs (m * m1)
#' @param threads number, the number of used threads for parallel process
#' @param lambda number, add a value on the diagonal of LD matrix
#' @param chisq number, generate a sparse LD matrix, if ld^2 * n < chisq, it will be set to 0
#' @param correlation logical, if TRUE, the LD matrix will be constructed by correlation of all pairs
#' @param out character, the path and prefix of output file, if not NULL, the LD matrix will stored in big.matrix, which will save much memory, but cost a bit more time
#' @param verbose logical, whether to print the log information

#' @examples
#' gwas_bfile_path <- system.file("extdata", "gwas_geno", package = "SumTool")
#' data <- read_plink(bfile=gwas_bfile_path, threads=1)
#' geno <- data$geno
#' ld <- LDcal(geno=geno, threads=1, verbose=FALSE)

LDcal <- function(geno = NULL, index = NULL, threads = 1, lambda = 0, chisq = 0, correlation = TRUE, out = NULL, verbose = TRUE)
{
	if(verbose)	version.info()
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Number of individuals:", nrow(geno), "\n")
	if(verbose)	cat("Number of SNPs:", ncol(geno), "\n")
	t1 <- as.numeric(Sys.time())
	if(is.null(out)){
		file.back <- FALSE
	}else{
		file.back <- TRUE
	}
	if(file.back){
		if(verbose)	cat("Transform genotype to matrix\n")
		geno <- as.matrix(geno)
		if(verbose)	cat("Write LD into bigmemory file\n")
		m <- ncol(geno)
		mc <- ifelse(is.null(index), m, length(index))
		txt1 <- ifelse(chisq > 0, ".sparse.bin", ".full.bin")
		txt2 <- ifelse(chisq > 0, ".sparse.desc", ".full.desc")
		backingfile <- paste(basename(out), txt1, sep="")
		descriptorfile <- paste(basename(out), txt2, sep="")
		if (file.exists(paste(dirname(out), "/", backingfile, sep=""))) file.remove(paste(dirname(out), "/", backingfile, sep=""))
		if (file.exists(paste(dirname(out), "/", descriptorfile, sep=""))) file.remove(paste(dirname(out), "/", descriptorfile, sep=""))
		bigmat <- bigmemory::filebacked.big.matrix(
			nrow = m,
			ncol = mc,
			type = "double",
			backingfile = backingfile,
			backingpath = dirname(out),
			descriptorfile = descriptorfile,
			dimnames = c(NULL, NULL)
		)
		SImpute_LD_bigm_c(geno, bigmat@address, index = index, chisq = chisq, lambda= lambda, haps = !correlation, threads = threads, verbose = verbose)
	}else{
		if(!is.big.matrix(geno))	stop("genotype should be in 'big.matrix' format!")
		bigmat <- SImpute_LD_norm_c(geno@address, index = index, chisq = chisq, lambda= lambda, haps = !correlation, threads = threads, verbose = verbose)
	}
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(bigmat)
}

#' Generalized Linear Model
#'
#' To do association using Generalized Linear Model 
#'
#' @param y vector, the phenotype vector with length n
#' @param geno big.matrix (n * m), genotype
#' @param map data.frame, the genomic information of SNPs. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2")
#' @param X matrix, covariates
#' @param threads number, the number of used threads for parallel process
#' @param verbose logical, whether to print the log information

#' @examples
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' geno <- data$geno
#' map <- data$map
#' y <- data$pheno[,1]
#' gwas <- LMreg(y=y, geno=geno, map=map, threads=1, verbose=FALSE)

LMreg <- function(y, geno = NULL, map = NULL, X = NULL, threads = 1, verbose = TRUE)
{
	if(verbose)	version.info()
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	if(is.null(geno))	stop("Please provide genotype!")
	if(is.null(map))	stop("Please provide map!")
	if(!is.big.matrix(geno))	stop("genotype should be in 'big.matrix' format!")
	if(ncol(geno) != nrow(map))	stop("Number of SNPs not equals between genotype and map!")
	if(ncol(map) != 5)	stop("Only 5 columns limited for map! (snp, chr, pos, a1, a2)")

	if(verbose)	cat("Number of individuals:", nrow(geno), "\n")
	if(is.matrix(y)){
		if(ncol(y) != 1)	stop("More than 1 phenotype!")
		y <- as.numeric(y)
	}
	if(nrow(geno) != length(y))	stop("Number of individuals not equals between phenotype and genotype!")
	if(verbose){
		if(sum(is.na(y)) != 0)	cat("Number of individuals with observations:", sum(!is.na(y)), "\n")
	}
	if(verbose)	cat("Number of SNPs:", ncol(geno), "\n")
	indx <- which(!is.na(y))
	yNa <- y[indx]
	if(!is.null(X)){
		X <- as.matrix(X)
		if(nrow(X) != length(y))	stop("Number of individuals not equals between phenotype and covariates!")
		if(verbose)	cat("Number of covariates:", ncol(X), "\n")
		X_col <- apply(X, 2, function(x) length(table(x)) > 1)
        X <- X[indx, X_col, drop=FALSE]
		X <- cbind(1, X)
	}else{
		X <- matrix(1, length(yNa))
	}
	if(verbose)	cat("Genome scanning...", "\n")
	t1 <- as.numeric(Sys.time())
	lmreg <- glm_c(y = yNa, X = X, indx = indx, geno@address, verbose = verbose, threads = threads)
	lmreg <- cbind(map, lmreg, length(yNa))
	colnames(lmreg) <- c("SNP", "Chr", "Pos", "A1", "A2", "Maf", "Beta", "SE", "P", "N")
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(lmreg)
}

#' LD score calculation
#'
#' To estimate LD score for each SNPs
#'
#' @param geno bigmemory (n * m), genotype coded as 0, 1, 2
#' @param map data.frame, the genomic information of SNPs. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2")
#' @param w int, size of windows in bp. Default is 1e6
#' @param b number, set the buffer size in bp of each window, default 125000
#' @param r2 logical, calculate r2 or r
#' @param adjust logical, whether to adjust the ldscore
#' @param verbose logical, whether to print the log information
#' @param threads number, the number of used threads for parallel process

#' @examples
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' 
#' # load data
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' geno <- data$geno
#' map <- data$map
#' ldscore <- LDscore(geno = geno, map = map, w = 100000, b=12500, threads = 1)

LDscore <- function(geno = NULL, map = NULL, w = 1000000, b = 500000, r2 = TRUE, adjust = TRUE, verbose = TRUE, threads = 1)
{

	if(verbose)	version.info()
	t1 <- as.numeric(Sys.time())
	
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	#check parameter
	if(verbose)	cat("Data and parameters check...")
	if(is.null(geno))	stop("Please provide geno!")
	if(!is.big.matrix(geno))	stop("genotype should be in 'big.matrix' format!")
	if(is.null(map))	stop("Please provide map!")
	if(ncol(map) != 5)	stop("Only 5 columns limited for map!")
	if(ncol(geno) != nrow(map))	stop("Number of SNPs not equals between geno and map!")
	if(!is.numeric(map[, 3]))	stop("Physical position in map should be numeric defined!")
	for(i in 1 : ncol(map)){
		if(is.factor(map[, i]))	map[, i] <- as.character.factor(map[, i])
	}
	#if(hasNA(geno@address, threads = threads))	stop("NA is not allowed in geno!")
	if(verbose) cat("(Qualified)\n")

	SNP_NA <- is.na(map[, 2 : 3])
	if(sum(SNP_NA) != 0){
		stop("NAs are not allowed in map!")
	}

	chr <- unique(map[, 2])

	#confirm parameters
	if(verbose)	cat("Number of chromosome is ", length(chr) , "\n", sep="")
	if(verbose)	cat("Number of total SNPs is ", ncol(geno) , "\n", sep="")
	if(verbose)	cat("Number of individuals is ", nrow(geno) , "\n", sep="")
	if(verbose)	cat("Window size is ", w / 1e6, "Mb", "\n", sep="")
	if(verbose)	cat("Buffer size is ", b / 1e3, "Kb", "\n", sep="")

	res <- NULL
	for(chri in chr){
		chri_index <- map[, 2] == chri
		if(verbose)	cat("Loop on chromosome", chri, "with", sum(chri_index), "SNPs\n")
		# chri_pos_min <- min(map[chri_index, 3])
		chri_pos_min <- 1
		chri_pos_max <- max(map[chri_index, 3])
		loop <- TRUE
		wind_min <- chri_pos_min
		wind_max <- wind_min + w
		while(loop){
			wind_min <- c(wind_min, wind_max[length(wind_max)] + 1)
			wind_max <- c(wind_max, wind_max[length(wind_max)] + w)
			if((wind_max[length(wind_max)] + b / 2) >= chri_pos_max)	loop <- FALSE
		}
		if(w >= chri_pos_max){
			wind_min <- wind_min[1]
			wind_max <- wind_max[1]
		}
		if(verbose)	cat("Total Number of windows: ", length(wind_min), "\n", sep="")
		if(verbose)	cat("Calculating LD score for windows...\n")
		wind_max[length(wind_max)] <- chri_pos_max
		wind_min_b <- wind_min - b / 2
		wind_max_b <- wind_max + b / 2
		wind_min_b[1] <- chri_pos_min
		wind_max_b[length(wind_max_b)] <- chri_pos_max
		mcf <- function(i){
			index1 <- which(chri_index & map[, 3] >= wind_min_b[i] & map[, 3] <= wind_max_b[i])
			if(length(index1) != 0){
				if(verbose)	cat("The ", i, "th window: Start[", min(map[index1, 3], na.rm = TRUE),"] ~ End[", max(map[index1, 3], na.rm = TRUE),"]\r", sep="")
				lds <- LDscore_c(geno@address, index1, r2 = r2, adjust = adjust, threads = threads, verbose = verbose)
				bin <- cbind(map[index1, ], lds)
				index <- (bin[, 3] >= wind_min[i]) & (bin[, 3] <= wind_max[i])
				bin <- bin[index, ]
			}else{
				bin <- NULL
			}
			return(bin)
		}
		chr_res <- lapply(1 : length(wind_min_b), mcf); gc()
		if(verbose)	cat("\n");
		res <- rbind(res, do.call(rbind, chr_res))
	}
	colnames(res) <- c("SNP", "Chr", "Pos", "A1", "A2", "Maf", "snp_num", "mean_rsq", "ldscore")
	res <- res[match(map[, 1], res[, 1]), ]
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Summary of LD Scores:\n")
	if(verbose)	print(summary(res$ldscore))
	if(verbose)	cat("MAF/LD Score Correlation Matrix:\n")
	if(verbose)	print(cor(res[, c("Maf", "ldscore")]))
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(res)
}

#' LD pruning
#'
#' To remove high correlated SNPs on base of r2 and MAF
#'
#' @param geno bigmemory (n * m), genotype coded as 0, 1, 2
#' @param map data.frame, the genomic information of SNPs. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2")
#' @param w int, size of windows in bp. Default is 1e6
#' @param b int, set the buffer size in bp of each window, default 125000
#' @param r2.cutoff double, the threshold of r2, smaller cutoff results in less remaining SNPs, default 0.25
#' @param verbose logical, whether to print the log information
#' @param threads int, the number of used threads for parallel process

#' @examples
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' 
#' # load data
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' geno <- data$geno
#' map <- data$map
#' snp <- LDprune(geno = geno, map = map, w = 100000, b=50000, threads = 1)

LDprune <- function(geno = NULL, map = NULL, w = 1000000, b = 500000, r2.cutoff = 0.25, verbose = TRUE, threads = 1)
{

	if(verbose)	version.info()
	t1 <- as.numeric(Sys.time())
	
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	#check parameter
	if(verbose)	cat("Data and parameters check...")
	if(is.null(geno))	stop("Please provide geno!")
	if(!is.big.matrix(geno))	stop("genotype should be in 'big.matrix' format!")
	if(is.null(map))	stop("Please provide map!")
	if(ncol(map) != 5)	stop("Only 5 columns limited for map!")
	if(ncol(geno) != nrow(map))	stop("Number of SNPs not equals between geno and map!")
	if(!is.numeric(map[, 3]))	stop("Physical position in map should be numeric defined!")
	for(i in 1 : ncol(map)){
		if(is.factor(map[, i]))	map[, i] <- as.character.factor(map[, i])
	}
	#if(hasNA(geno@address, threads = threads))	stop("NA is not allowed in geno!")
	if(verbose) cat("(Qualified)\n")

	SNP_NA <- is.na(map[, 2 : 3])
	if(sum(SNP_NA) != 0){
		stop("NAs are not allowed in map!")
	}

	chr <- unique(map[, 2])

	#confirm parameters
	if(verbose)	cat("Number of chromosome is ", length(chr) , "\n", sep="")
	if(verbose)	cat("Number of total SNPs is ", ncol(geno) , "\n", sep="")
	if(verbose)	cat("Number of individuals is ", nrow(geno) , "\n", sep="")
	if(verbose)	cat("Window size is ", w / 1e6, "Mb", "\n", sep="")
	if(verbose)	cat("Buffer size is ", b / 1e3, "Kb", "\n", sep="")
	if(verbose)	cat("r2 threshold ", r2.cutoff, "\n", sep="")

	res <- NULL
	for(chri in chr){
		chri_index <- map[, 2] == chri
		if(verbose)	cat("Loop on chromosome", chri, "with", sum(chri_index), "SNPs\n")
		# chri_pos_min <- min(map[chri_index, 3])
		chri_pos_min <- 1
		chri_pos_max <- max(map[chri_index, 3])
		loop <- TRUE
		wind_min <- chri_pos_min
		wind_max <- wind_min + w
		while(loop){
			wind_min <- c(wind_min, wind_max[length(wind_max)] + 1)
			wind_max <- c(wind_max, wind_max[length(wind_max)] + w)
			if((wind_max[length(wind_max)] + b / 2) >= chri_pos_max)	loop <- FALSE
		}
		if(w >= chri_pos_max){
			wind_min <- wind_min[1]
			wind_max <- wind_max[1]
		}
		if(verbose)	cat("Total Number of windows: ", length(wind_min), "\n", sep="")
		wind_max[length(wind_max)] <- chri_pos_max
		wind_min_b <- wind_min - b / 2
		wind_max_b <- wind_max + b / 2
		wind_min_b[1] <- chri_pos_min
		wind_max_b[length(wind_max_b)] <- chri_pos_max
		# mcf <- function(i){
		chr_res <- NULL
		for(i in 1 : length(wind_min)){
			index1 <- which(chri_index & map[, 3] >= wind_min_b[i] & map[, 3] <= wind_max_b[i])
			if(length(index1) != 0){
				if(verbose)	cat("The ", i, "th window: Start[", min(map[index1, 3], na.rm = TRUE),"] ~ End[", max(map[index1, 3], na.rm = TRUE),"]\r", sep="")
				snp_indx <- LDprune_c(geno@address, index1, r2_cutoff = r2.cutoff, threads = threads, verbose = verbose)
				bin <- map[index1[snp_indx == 1], ]
				index <- (bin[, 3] >= wind_min[i]) & (bin[, 3] <= wind_max[i])
				# bin <- bin[index, 1]
				chr_res <- c(chr_res, bin[index, 1])
			}else{
				# bin <- NULL
			}
			# return(bin)
		}
		# chr_res <- unlist(lapply(1 : length(wind_min), mcf)); gc()
		if(verbose)	cat("\n");
		if(verbose)	cat("Pruned", sum(chri_index) - length(chr_res), "variants from chromosome", chri, "with", length(chr_res), "left\n")
		res <- c(res, chr_res)
	}
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("After pruned, total", length(res), "remains\n")
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(res)
}

#' LD clumping
#'
#' To remove high correlated SNPs on base of r2 and p-values
#'
#' @param geno bigmemory (n * m), genotype coded as 0, 1, 2
#' @param map data.frame, the genomic information of SNPs. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2")
#' @param p data.frame, at least 2 columns, the first column should be SNP names, the last column should be the pvalues
#' @param w int, size of windows in bp. Default is 1e6
#' @param r2.cutoff double, the threshold of r2, smaller cutoff results in less remaining SNPs, default 0.25
#' @param p.cutoff double, the threshold of p-values, smaller threshold results in less remaining SNPs, default 1
#' @param verbose logical, whether to print the log information
#' @param threads int, the number of used threads for parallel process

#' @examples
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' p_path <- system.file("extdata", "P.txt", package = "SumTool")

#' # load data
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' geno <- data$geno
#' map <- data$map
#' pdata <- read.table(p_path, header = TRUE)
#' snp <- LDclump(geno=geno, map=map, p=pdata, p.cutoff=1, r2.cutoff=0.25, w=100000)

LDclump <- function(geno = NULL, map = NULL, p = NULL, w = 1000000, r2.cutoff = 0.25, p.cutoff = 1, verbose = TRUE, threads = 1)
{

	if(verbose)	version.info()
	t1 <- as.numeric(Sys.time())
	
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	#check parameter
	if(verbose)	cat("Data and parameters check...")
	if(is.null(geno))	stop("Please provide geno!")
	if(!is.big.matrix(geno))	stop("genotype should be in 'big.matrix' format!")
	if(is.null(map))	stop("Please provide map!")
	if(is.null(p))	stop("Please provide summary statistic with p-values included at the end of column!")
	if(ncol(p) < 2)	stop("At least 2 columns should be provided, the first column should be SNP names, the last column should be the pvalues")
	if(!is.numeric(p[, ncol(p)]))	stop("p-values should be numeric defined!")
	if(ncol(map) != 5)	stop("Only 5 columns limited for map!")
	if(ncol(geno) != nrow(map))	stop("Number of SNPs not equals between geno and map!")
	if(!is.numeric(map[, 3]))	stop("Physical position in map should be numeric defined!")
	for(i in 1 : ncol(map)){
		if(is.factor(map[, i]))	map[, i] <- as.character.factor(map[, i])
	}
	#if(hasNA(geno@address, threads = threads))	stop("NA is not allowed in geno!")
	if(verbose) cat("(Qualified)\n")

	SNP_NA <- is.na(map[, 2 : 3])
	if(sum(SNP_NA) != 0){
		stop("NAs are not allowed in map!")
	}

	mergedSNP <- intersect(map[, 1], p[, 1])
	if(length(mergedSNP) == 0)	stop("No shared SNPs between 'p' and 'map'!")
	
	#confirm parameters
	if(verbose)	cat("Number of chromosome in geno ", length(unique(map[, 2])), "\n", sep="")
	if(verbose)	cat("Number of total SNPs in geno ", ncol(geno), "\n", sep="")
	if(verbose)	cat("Number of individuals in geno ", nrow(geno), "\n", sep="")
	if(verbose)	cat("Number of SNPs with p-values ", nrow(p), "\n", sep="")
	if(verbose)	cat("Number of shared SNPs ", length(mergedSNP), "\n", sep="")
	if(verbose)	cat("Window size is ", w / 1e6, "Mb", "\n", sep="")
	if(verbose)	cat("r2 threshold ", r2.cutoff, "\n", sep="")
	if(verbose)	cat("P-value threshold ", p.cutoff, "\n", sep="")
	p <- p[match(mergedSNP, p[, 1]), ]

	# remove snps on base of pvalue
	p <- p[p[, ncol(p)] <= p.cutoff, ]
	if(verbose)	cat("Remove ", length(mergedSNP) - nrow(p), " SNPs on base of p-value threshold\n", sep="")
	if(nrow(p) == 0)	stop("No SNPs left, please reset 'p.cutoff'!")

	p.map <- map[match(p[, 1], map[, 1]), ]
	chr <- unique(p.map[, 2])

	res <- NULL
	for(chri in chr){
		chri_index <- p.map[, 2] == chri
		if(verbose)	cat("Loop on chromosome", chri, "with", sum(chri_index), "SNPs\n")
		# chri_pos_min <- min(p.map[chri_index, 3])
		chri_pos_min <- 1
		chri_pos_max <- max(p.map[chri_index, 3])
		loop <- TRUE
		wind_min <- chri_pos_min
		wind_max <- wind_min + w
		while(loop){
			wind_min <- c(wind_min, wind_max[length(wind_max)] + 1)
			wind_max <- c(wind_max, wind_max[length(wind_max)] + w)
			if((wind_max[length(wind_max)]) >= chri_pos_max)	loop <- FALSE
		}
		if(w >= chri_pos_max){
			wind_min <- wind_min[1]
			wind_max <- wind_max[1]
		}
		if(verbose)	cat("Total Number of windows: ", length(wind_min), "\n", sep="")
		wind_max[length(wind_max)] <- chri_pos_max
		mcf <- function(i){
			index1 <- which(chri_index & p.map[, 3] >= wind_min[i] & p.map[, 3] <= wind_max[i])
			if(length(index1) != 0){
				if(verbose)	cat("The ", i, "th window: Start[", min(p.map[index1, 3], na.rm = TRUE),"] ~ End[", max(p.map[index1, 3], na.rm = TRUE),"]\r", sep="")
				snp_indx <- LDclump_c(geno@address, match(p.map[index1, 1], map[, 1]), r2_cutoff = r2.cutoff, p = p[index1, ncol(p)], threads = threads, verbose = verbose)
				bin <- p.map[index1[snp_indx == 1], 1]
			}else{
				bin <- NULL
			}
			return(bin)
		}
		chr_res <- unlist(lapply(1 : length(wind_min), mcf)); gc()
		if(verbose)	cat("\n");
		if(verbose)	cat("Clumped", sum(chri_index) - length(chr_res), "variants from chromosome", chri, "with", length(chr_res), "left\n")
		res <- c(res, chr_res)
	}
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("After clumped, total", length(res), "remains\n")
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(res)
}

#' Summary data based SNP BLUP
#'
#' To estimate joint SNP effect using GWAS summary statistics
#'
#' @param sumstat data.frame, GWAS summary statistics of a trait. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2", "BETA", "SE, "N")
#' @param geno bigmemory (n * m), genotype coded as 0, 1, 2
#' @param map data.frame, the genomic information of SNPs. The columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2")
#' @param lambda double, the ridge regression coefficient, lambda=m*(1/h2-1), m is the number of SNPs, h2 is the heritability of trait
#' @param w int, size of windows in bp. Default is 1e6
#' @param threads number, the number of used threads for parallel process
#' @param verbose logical, whether to print the log information

#' @examples
#' sumstat_path <- system.file("extdata", "typed.marginal", package = "SumTool")
#' ref_bfile_path <- system.file("extdata", "ref_geno", package = "SumTool")
#' 
#' # load data
#' sumstat <- read.table(sumstat_path, header=TRUE)
#' data <- read_plink(bfile=ref_bfile_path, threads=1)
#' geno <- data$geno
#' map <- data$map
#' h2 <- 0.5
#' lambda = nrow(sumstat)*(1/h2-1)
#' eff <- SBLUP(sumstat=sumstat, geno=geno, map=map, lambda=lambda, threads=1)

SBLUP <- function(sumstat = NULL, geno = NULL, map = NULL, lambda = NULL, w = 1e6, threads = 1, verbose = TRUE){

	if(verbose)	version.info()
	t1 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	#check parameter
	if(verbose)	cat("Data and parameters check...")
	if(is.null(sumstat))	stop("Please provide 'sumstat'!")
	if(!is.big.matrix(geno))	stop("genotype should be in 'big.matrix' format!")
	if(is.null(geno))	stop("Please provide 'geno'!")
	if(is.null(map))	stop("Please provide 'map'!")
	if(ncol(map) != 5)	stop("Only 5 columns limited for map! (SNP Chr Pos A1 A2)")
	if(ncol(sumstat) != 8)	stop("Only 8 columns limited for typed! (SNP Chr Pos A1 A2 beta se N)")
	if(!is.numeric(map[, 3]))	stop("Physical position in map should be numeric defined!")
	if(!is.numeric(sumstat[, 3]))	stop("Physical position in sumstat should be numeric defined!")
	for(i in 1 : ncol(map)){
		if(is.factor(map[, i]))	map[, i] <- as.character.factor(map[, i])
		if(is.factor(sumstat[, i]))	sumstat[, i] <- as.character.factor(sumstat[, i])
	}
	if(sum(is.na(map[, 3])) != 0)	stop(paste("NA is not allowed for physical position of map!", sep=""))
	if(sum(is.na(sumstat[, 3])) != 0)	stop(paste("NA is not allowed for physical position of sumstat!", sep=""))
	if(sum(is.na(sumstat[, -c(1:5)])) != 0)	stop(paste("NA is not allowed in sumstat information!", sep=""))
	if(ncol(geno) != nrow(map))	stop("Number of SNPs not equals between geno and map!")
	if(any(duplicated(map[, 1])))	stop("duplicated SNP names exist in map file!")
	if(any(duplicated(sumstat[, 1])))	stop("duplicated SNP names exist in sumstat file!")
	#if(hasNA(geno@address, threads = threads))	stop("NA is not allowed in geno!")
	if(is.null(lambda))	stop("'lambda should be provided! \nlambda=m*((1/h2)-1), where m is the number of SNPs, h2 is the heritability of trait.")
	if(verbose) cat("(Qualified)\n")

	#confirm parameters
	if(verbose)	cat("Number of total SNPs in geno is ", ncol(geno) , "\n", sep="")
	if(verbose)	cat("Number of individuals in geno is ", nrow(geno) , "\n", sep="")
	if(verbose)	cat("Number of total SNPs in summary statistics is ", nrow(sumstat) , "\n", sep="")
	N <- round(mean(sumstat[, ncol(sumstat)]))
	if(verbose)	cat("Number of individuals in summary statistics is ", N, "\n", sep="")
	mergedSNP <- intersect(sumstat[, 1], map[, 1])
	if(length(mergedSNP) == 0)	stop("No shared SNPs between sumstat and map!")
	sumstat <- sumstat[match(mergedSNP, sumstat[, 1]), ]
	row_logi <- NULL
	for(i in 2 : 5){
		row_logi <- cbind(row_logi, sumstat[, i] == map[match(mergedSNP, map[, 1]), i])
	}
	indx <- apply(row_logi, 1, all); rm(row_logi); gc()
	if(!all(indx)){
		sumstat <- sumstat[indx, ]
		if(verbose)	cat("(Warning:", sum(!indx), "SNPs sharing same names between reference and sumstat don't have equal map information!)\n")
	}
	if(nrow(sumstat) == 0)	stop("No shared SNPs between sumstat and map! Please check the chr, pos, A1 and A2 for the SNPs with same names")
	if(verbose)	cat("Ridge regression coefficient", lambda, "\n")
	if(verbose)	cat("Number of shared typed SNPs is ", nrow(sumstat), "\n", sep="")
	if(verbose)	cat("Window size is ", w / 1e6, "Mb", "\n", sep="")
	if(verbose)	cat("Estimating Joint effect...\n") 
	chr <- unique(sumstat[, 2])

	res <- NULL
	for(chri in chr){
		chri_index <- sumstat[, 2] == chri
		if(verbose)	cat("Loop on chromosome", chri, "with", sum(chri_index), "SNPs\n")
		# chri_pos_min <- min(sumstat[chri_index, 3])
		chri_pos_min <- 1
		chri_pos_max <- max(sumstat[chri_index, 3])
		loop <- TRUE
		wind_min <- chri_pos_min
		wind_max <- chri_pos_min + w
		while(loop){
			wind_min <- c(wind_min, wind_max[length(wind_max)] + 1)
			wind_max <- c(wind_max, wind_max[length(wind_max)] + w)
			if((wind_max[length(wind_max)]) >= chri_pos_max)	loop <- FALSE
		}
		if(w >= chri_pos_max){
			wind_min <- wind_min[1]
			wind_max <- wind_max[1]
		}
		wind_max[length(wind_max)] <- chri_pos_max
		if(verbose)	cat("Total Number of windows: ", length(wind_min), "\n", sep="")
		mcf <- function(i){
			index1 <- which(chri_index & sumstat[, 3] >= wind_min[i] & sumstat[, 3] <= wind_max[i])
			index2 <- match(sumstat[index1, 1], map[, 1])
			if(length(index1) != 0){
				if(verbose)	cat("The ", i, "th window: Start[", min(sumstat[index1, 3], na.rm = TRUE),"] ~ End[", max(sumstat[index1, 3], na.rm = TRUE),"]\r", sep="")
				effect <- sblup_bin(geno@address, N, index2, sumstat[index1, 6], lambda = lambda, verbose = verbose, threads = threads)
				effect <- as.vector(effect)
			}else{
				effect <- NULL
			}
			return(effect)
		}
		chr_res <- lapply(1 : length(wind_min), mcf); gc()
		if(verbose)	cat("\n");
		res <- c(res, unlist(chr_res))
	}
	res <- cbind(sumstat[, c(1 : 5)], res)
	colnames(res) <- c("SNP", "Chr", "Pos", "A1", "A2", "Effect")
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(res)
}

#' LD regression
#'
#' To estimate heritability or genetic correlation
#'
#' @param sumstat data.frame or list, results of single trait or multiple traits of GWAS summary statistics. For each summary statistic of a trait, the columns should be in the order of c("SNP", "Chr", "Pos", "A1", "A2", "BETA", "SE", "N"). If it's a data.frame, the heritability will be estimated, if it's a list containing multiple GWAS summary statistics, both heritability of traits and the genetic correlation of pair of traits will be estimated.
#' @param ldscore data.frame, the ldscore for each SNP, the columns should be prepared in order of c("SNP", "Chr", "Pos", "A1", "A2", "MAF", "ldscore")
#' @param wld data.frame, at least 2 colmuns, the first column is th SNP names, the last column is the weight for each SNP. Default equals to the 'ldscore' option
#' @param maxz2 double, max chi^2. Default is 80
#' @param maf double, Minor allele frequency lower bound. Default is MAF > 0.05
#' @param nblock int, number of blocks for jackknife estimation
#' @param verbose logical, whether to print the log information

#' @examples
#' # single trait 
#' sumstat1_path <- system.file("extdata", "sumstat1", package = "SumTool")
#' ldscore_path <- system.file("extdata", "ldscore", package = "SumTool")
#' sumstat1 <- read.table(sumstat1_path, header=TRUE)
#' ldscore <- read.table(ldscore_path, header=TRUE)
#' res1 <- LDreg(sumstat = sumstat1, ldscore = ldscore)
#' 
#' #multiple traits
#' sumstat2_path <- system.file("extdata", "sumstat2", package = "SumTool")
#' sumstat2 <- read.table(sumstat2_path, header=TRUE)
#' res2 <- LDreg(sumstat = list(sumstat1, sumstat2), ldscore = ldscore)

LDreg <- function(sumstat = NULL, ldscore = NULL, wld = NULL, maxz2 = 80, maf = 0.05, nblock = 200, verbose = TRUE){

	if(verbose)	version.info()
	if(verbose)	cat("Analysis started:", as.character(Sys.time()), "\n")
	t1 <- as.numeric(Sys.time())

	h2_est <- function(sumstat, ldscore, wld = NULL, maxz2, maf, nblock, verbose){
		if(ncol(sumstat) != 8)	stop("Number of columns for sumstat should be 8 (snp chr pos a1 a2 beta se n).")
		if(ncol(ldscore) != 7)	stop("Number of columns for ldscore should be 7 (snp chr pos a1 a2 maf ldscore).")
		if(verbose)	cat("Number of SNPs for summary statistics", nrow(sumstat), "\n")
		if(verbose)	cat("Number of SNPs for ldscore", nrow(ldscore), "\n")
		if(sum(is.na(ldscore[, ncol(ldscore)])) != 0)	stop("NA is not allowed!")
		if(sum(is.na(sumstat[, -c(1:5)])) != 0)	stop("NA is not allowed!")
		if(any(duplicated(ldscore[, 1])))	stop("duplicated SNP names exist in ldscore file!")
		if(any(duplicated(sumstat[, 1])))	stop("duplicated SNP names exist in sumstat file!")

		# order ldscore by the positin on genome
		ldscore <- ldscore[order(ldscore[, 2], ldscore[, 3]), ]
		M <- sum(ldscore[, 6] > maf)
		
		if(!is.null(wld)){
			if(ncol(as.matrix(wld)) < 2)	stop("At least 2 columns should be provided!")
			if(sum(is.na(wld[, ncol(wld)])) != 0)	stop("NA is not allowed in wld!")
			if(any(duplicated(wld[, 1])))	stop("duplicated SNP names exist in wld!")
			if(nrow(ldscore) != nrow(wld))	stop("Number of SNPs between ldscore and wld not match!")
			mergedSNP <- intersect(ldscore[, 1], wld[, 1])
			if(length(mergedSNP) == 0)	stop("No shared SNPs between ldscore and wld!")
			if(verbose)	cat("After merging with regression weight", length(mergedSNP), "SNPs remain", "\n")
			if(length(mergedSNP) == nrow(ldscore)){
				wld <- wld[match(ldscore[, 1], wld[, 1]), ]
			}else{
				ldscore <- ldscore[ldscore[, 1] %in% mergedSNP, ]
				wld <- wld[match(ldscore[, 1], wld[, 1]), ]
			}
		}else{
			wld <- ldscore
		}

		#genomic inflation factor
		sumstat <- cbind(sumstat, (sumstat[, 6] / sumstat[, 7])^2)
		lambda <- median(sumstat[, ncol(sumstat)]) / qchisq(1/2, df = 1,lower.tail=F)

		mergedSNP <- intersect(ldscore[, 1], sumstat[, 1])
		if(length(mergedSNP) == 0)	stop("No shared SNPs between ldscore and sumstat!")
		if(length(mergedSNP) != nrow(ldscore)){
			indx <- ldscore[, 1] %in% mergedSNP
			ldscore <- ldscore[indx, ]
			wld <- wld[indx, ]
		}
		sumstat <- sumstat[match(ldscore[, 1], sumstat[, 1]), ]
		row_logi <- NULL
		for(i in 2 : 5){
			row_logi <- cbind(row_logi, ldscore[, i] == sumstat[, i])
		}
		indx <- apply(row_logi, 1, all); rm(row_logi); gc()
		if(!all(indx)){
			if(verbose)	cat("(Warning:", sum(!indx), "SNPs sharing same names between sumstat and ldscore don't have equal map information!)\n")
			ldscore <- ldscore[indx, ]
			wld <- wld[indx, ]
			sumstat <- sumstat[indx, ]
		}
		if(verbose)	cat("After merging with reference panel", nrow(sumstat), "SNPs remain", "\n")
		N <- mean(sumstat[, 8])
		if(verbose)	cat("Total", M, "SNPs that MAF >", maf, "\n")
		if(verbose)	cat("Total", N, "individuals included", "\n")
		if(verbose)	cat("Using two-step estimator with cutoff at", maxz2, "\n")
		if(nblock > nrow(sumstat))	nblock <- nrow(sumstat)
		h2 <- ldreg_h2(z2 = sumstat[, ncol(sumstat)], ld = ldscore[, ncol(ldscore)], wld = wld[, ncol(wld)], M = M, N = N, maxz2 = maxz2, nblock = nblock)
		h2 <- as.vector(h2)
		if(verbose)	cat("Estimated h2: ", h2[1], " (",h2[2], ")", sep="", "\n")
		if(verbose)	cat("Lambda GC:", lambda, "\n")
		if(verbose)	cat("Mean Chi^2:", mean(sumstat[, ncol(sumstat)]), "\n")
		if(verbose)	cat("Intercept: ", h2[3], " (",h2[4], ")", sep="", "\n")
		names(h2) <- c("Heritability", "SE", "Intercept", "SE")
		return(h2)
	}

	l <- ldscore[, ncol(ldscore)]
	l[l < 1] <- 1
	ldscore[, ncol(ldscore)] <- l
	if(!is.null(wld)){
		l <- wld[, ncol(wld)]
		l[l < 1] <- 1
		wld[, ncol(wld)] <- l	
	}
	if(is.matrix(sumstat) || is.data.frame(sumstat)){

		#heritability estimation
		est <- h2_est(sumstat = sumstat, ldscore = ldscore, wld = wld, maxz2 = maxz2, maf = maf, nblock = nblock, verbose = verbose)
	}else if(is.list(sumstat)){
		if(length(sumstat) == 1){
			
			#heritability estimation
			est <- h2_est(sumstat = sumstat, ldscore = ldscore, wld = wld, maxz2 = maxz2, maf = maf, nblock = nblock, verbose = verbose)
		}else{
			if(verbose)	cat("Bivariate LD regression on", length(sumstat), "traits\n")
			est <- NULL
			rg <- diag(1, length(sumstat))
			colnames(rg) <- paste0("trait", 1 : length(sumstat))
			rownames(rg) <- colnames(rg)

			#heritability estimation
			for(i in 1 : length(sumstat)){
				if(verbose)	cat("\nHeritability of phenotype", i, "\n")
				if(verbose)	cat("--------------------------------------------\n")
				est0 <- h2_est(sumstat = sumstat[[i]], ldscore = ldscore, wld = wld, maxz2 = maxz2, maf = maf, nblock = nblock, verbose = verbose)
				est <- rbind(est, est0)
				rownames(est)[i] <- paste0("h2_", i)
			}

			# order ldscore by the positin on genome
			ldscore <- ldscore[order(ldscore[, 2], ldscore[, 3]), ]
			M <- sum(ldscore[, 6] > maf)

			if(!is.null(wld)){
				mergedSNP <- intersect(ldscore[, 1], wld[, 1])
				if(length(mergedSNP) == nrow(ldscore)){
					wld <- wld[match(ldscore[, 1], wld[, 1]), ]
				}else{
					ldscore <- ldscore[ldscore[, 1] %in% mergedSNP, ]
					wld <- wld[match(ldscore[, 1], wld[, 1]), ]
				}
			}else{
				wld <- ldscore
			}
			ldscore_ori <- ldscore
			wld_ori <- wld

			#genetic correlation estimation
			for(i in 1 : (length(sumstat) - 1)){
				sumstat1 <- sumstat[[i]]
				ldscore <- ldscore_ori
				wld <- wld_ori

				mergedSNP <- intersect(ldscore[, 1], sumstat1[, 1])
				if(length(mergedSNP) != nrow(ldscore)){
					indx <- ldscore[, 1] %in% mergedSNP
					ldscore <- ldscore[indx, ]
					wld <- wld[indx, ]
				}
				sumstat1 <- sumstat1[match(ldscore[, 1], sumstat1[, 1]), ]
				row_logi <- NULL
				for(k in 2 : 5){
					row_logi <- cbind(row_logi, ldscore[, k] == sumstat1[, k])
				}
				indx <- apply(row_logi, 1, all); rm(row_logi); gc()
				if(!all(indx)){
					ldscore <- ldscore[indx, ]
					wld <- wld[indx, ]
					sumstat1 <- sumstat1[indx, ]
				}

				for(j in (i + 1) : length(sumstat)){
					sumstat2 <- sumstat[[j]]

					mergedSNP <- intersect(sumstat1[, 1], sumstat2[, 1])
					if(length(mergedSNP) == 0)	stop("No shared SNPs between sumstat", i, " and sumstat", j, "!")
					if(length(mergedSNP) != nrow(sumstat1)){
						indx <- sumstat1[, 1] %in% mergedSNP
						sumstat1 <- sumstat1[indx, ]
						ldscore <- ldscore[indx, ]
						wld <- wld[indx, ]
					}
					if(length(mergedSNP) != nrow(sumstat2)){
						sumstat2 <- sumstat2[sumstat2[, 1] %in% mergedSNP, ]
					}
					sumstat2 <- sumstat2[match(mergedSNP, sumstat2[, 1]), ]

					if(verbose)	cat("\nGenetic Covariance between phenotype", i, "and", j, "\n")
					if(verbose)	cat("--------------------------------------------\n")
					if(verbose)	cat("Number of shared SNPs:", length(mergedSNP), "\n")
					if(verbose)	cat("Number of individuals for two traits:", mean(sumstat1[, 8]), mean(sumstat2[, 8]), "\n")
					if(nblock > nrow(sumstat1)){
						nblock1 <- nrow(sumstat1)
					}else{
						nblock1 <- nblock
					}
					est0 <- ldreg_rg((sumstat1[, 6] / sumstat1[, 7]), (sumstat2[, 6] / sumstat2[, 7]), ld = ldscore[, ncol(ldscore)], wld = wld[, ncol(wld)], M = M, N1 = mean(sumstat1[, 8]), N2 = mean(sumstat2[, 8]), nblock1)
					est0 <- as.vector(est0)
					if(verbose)	cat("Estimated gencov: ", est0[1], " (",est0[2], ")", sep="", "\n")
					if(verbose)	cat("Mean z1*z2:", mean((sumstat1[, 6] / sumstat1[, 7]) * (sumstat2[, 6] / sumstat2[, 7])), "\n")
					if(verbose)	cat("Intercept: ", est0[3], " (",est0[4], ")", sep="", "\n")
					rg[i, j] <- est0[1] / sqrt(est[i, 1] * est[j, 1]) -> rg[j, i]
					est <- rbind(est, est0)
					rownames(est)[nrow(est)] <- paste0("h2_", i, "_", j)
				}
			}
			if(verbose)	cat("Genetic Correlation:\n")
			if(verbose)	print(rg)
		}
	}else{
		stop("Unknow type of format!")
	}
	t2 <- as.numeric(Sys.time())
	if(verbose)	cat("Analysis finished:", as.character(Sys.time()), "\n")
	if(verbose)	cat("Total Running time:", times(t2-t1), "\n")
	return(est)
}
