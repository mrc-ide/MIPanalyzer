
#------------------------------------------------
#' @useDynLib MIPanalyzer, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
NULL

#------------------------------------------------
#' @title Convert vcf to biallelic mipanalyzer data class
#'
#' @description Convert vcf to biallelic mipanalyzer data class.
#'
#' @param file path to vcf file.
#' @param vcfR object of class \code{vcfR}.
#' @param verbose if reading from file, whether to read in verbose manner.
#'
#' @export

vcf2mipanalyzer_biallelic <- function(file = NULL, vcfR = NULL, verbose = TRUE) {
  
  # check inputs
  if (!xor(!is.null(file), !is.null(vcfR))) {
    stop("Must specify one input: either a raw vcf file path or a vcfR object")
  }
  
  # get vcf object
  if(!is.null(vcfR)){
    assert_custom_class(vcfR, "vcfR")
    vcf <- vcfR
  } else {
    assert_file_exists(file)
    vcf <- vcfR::read.vcfR(file = file, verbose = verbose)
  }
  
  # check that vcf is normalised
  vcf_unnormalised <- vcf@fix %>%
                      tibble::as.tibble(.) %>%
                      dplyr::group_by(CHROM, POS) %>%
                      dplyr::summarise(locicount = n()) %>%
                      dplyr::ungroup() %>%
                      dplyr::summarise(maxlocicount = max(locicount)) != 1
  if (vcf_unnormalised) {
    stop("This is not a normalized vcf. Consider running bcftools norm, and/or review how the vcf was created")
  }
  
  # check that biallelic at all loci
  if (!all(vcfR::is.biallelic(vcf))) {
    stop("All loci must be biallelic")
  }
  
  # extract coverage and counts matrices
  coverage <- t(vcfR::extract.gt(vcf, element = "DP", as.numeric = T))
  counts_raw <- t(vcfR::extract.gt(vcf, element = "AD"))
  counts <- masplit(counts_raw, record = 1, sort = FALSE, decreasing = FALSE)
  
  # check that all missing fields correspond between coverage and counts
  if (!identical(is.na(coverage), is.na(counts))) {
    stop("Cannot have fields that are NA in coverage but not in counts, or vice versa")
  }
  
  # extract sample info
  samples <- data.frame(SAMPLE_ID = colnames(vcf@gt)[2:ncol(vcf@gt)], stringsAsFactors = FALSE)
  
  # extract loci and specify some columns classes
  loci <- as.data.frame(vcf@fix)
  loci$POS <- as.numeric(as.character(loci$POS))
  loci$QUAL <- as.numeric(as.character(loci$QUAL))
  
  # initialise filter history
  filter_history <- data.frame(description = "raw data",
                               samples = nrow(coverage),
                               loci = ncol(coverage),
                               n_missing = sum(is.na(coverage)),
                               prop_missing = mean(is.na(coverage)),
                               function_call = NA,
                               stringsAsFactors = FALSE)
  
  # create return list
  ret <- list(coverage = coverage,
              counts = counts,
              samples = samples,
              loci = loci,
              filter_history = filter_history)
  
  # return in mipanalyzer_biallelic class
  class(ret) <- "mipanalyzer_biallelic"
  return(ret)
}

#------------------------------------------------
#' @title Filter out some samples
#'
#' @description Filter out some samples.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param sample_filter boolean vector specifying whether to keep (\code{TRUE})
#'   or drop (\code{FALSE}) each sample.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_samples <- function(x, sample_filter, description = "") {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_vector(sample_filter)
  assert_logical(sample_filter)
  assert_eq(length(sample_filter), nrow(x$samples))
  
  # apply filter
  x$coverage <- x$coverage[sample_filter,]
  x$counts <- x$counts[sample_filter,]
  x$samples <- x$samples[sample_filter,]
  
  # record filter
  function_call <- paste0(deparse(match.call()), collapse = "")
  x$filter_history <- rbind(x$filter_history, list(description = description,
                                                   samples = nrow(x$coverage),
                                                   loci = ncol(x$coverage),
                                                   n_missing = sum(is.na(x$coverage)),
                                                   prop_missing = mean(is.na(x$coverage)),
                                                   function_call = function_call))
  
  return(x)
}

#------------------------------------------------
#' @title Filter out some loci
#'
#' @description Filter out some loci.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param locus_filter boolean vector specifying whether to keep (\code{TRUE})
#'   or drop (\code{FALSE}) each locus.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_loci <- function(x, locus_filter, description = "") {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_vector(locus_filter)
  assert_logical(locus_filter)
  assert_eq(length(locus_filter), nrow(x$loci))
  
  # apply filter
  x$coverage <- x$coverage[,locus_filter]
  x$counts <- x$counts[,locus_filter]
  x$loci <- x$loci[locus_filter,]
  
  # record filter
  function_call <- paste0(deparse(match.call()), collapse = "")
  x$filter_history <- rbind(x$filter_history, list(description = description,
                                                   samples = nrow(x$coverage),
                                                   loci = ncol(x$coverage),
                                                   n_missing = sum(is.na(x$coverage)),
                                                   prop_missing = mean(is.na(x$coverage)),
                                                   function_call = function_call))
  
  return(x)
}

#------------------------------------------------
#' @title Filter out over-counts
#'
#' @description Filter out over-counts, defined as count > coverage.
#'   Replace any such element with NA.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_overcounts <- function(x, description = "replace overcounts with NA") {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  
  # replace over-counts with NA
  w <- which(x$counts > x$coverage, arr.ind = TRUE)
  x$coverage[w] <- NA
  x$counts[w] <- NA
  
  # record filter
  function_call <- paste0(deparse(match.call()), collapse = "")
  x$filter_history <- rbind(x$filter_history, list(description = description,
                                                   samples = nrow(x$coverage),
                                                   loci = ncol(x$coverage),
                                                   n_missing = sum(is.na(x$coverage)),
                                                   prop_missing = mean(is.na(x$coverage)),
                                                   function_call = function_call))
  
  return(x)
}

#------------------------------------------------
#' @title Explore sample coverage prior to filtering
#'
#' @description See \code{filter_coverage_samples()} function first. Explore
#'   values of \code{min_coverage} and \code{max_missing} without applying any
#'   filtering.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param min_coverage any coverage value \code{< min_coverage} is replaced by \code{NA}.
#' @param max_missing any sample with \code{> max_missing} values equal \code{NA} are dropped.
#' @param breaks number of breaks spanning the range \code{[0,100]}.
#'
#' @export

explore_filter_coverage_samples <- function(x,
                                            min_coverage = 5,
                                            max_missing = 50,
                                            breaks = 100) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_missing)
  assert_bounded(max_missing, right = 100)
  
  # get percent missing per sample
  x$coverage[x$coverage < min_coverage] <- NA
  percent_missing <- rowMeans(is.na(x$coverage)) * 100
  percent_drop <- round(mean(percent_missing > max_missing)*100, digits = 2)
  
  # construct main title
  main_title <- paste0("min_coverage = ", min_coverage,
                       "\nmax_missing = ", max_missing,
                       "\nsamples dropped = ", percent_drop, "%")
  
  # produce plot
  plot1 <- ggplot2::ggplot(data = data.frame(x = percent_missing)) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_histogram(ggplot2::aes(x = x), fill = "#4575B4", breaks = seq(0,100,l = breaks))
  plot1 <- plot1 + ggplot2::geom_vline(xintercept = max_missing, linetype = "dashed")
  plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(0,100))
  plot1 <- plot1 + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0))
  plot1 <- plot1 + ggplot2::xlab("% missing loci per sample")
  plot1 <- plot1 + ggplot2::ggtitle(main_title)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Filter samples based on coverage
#'
#' @description Set a coverage threshold: any coverage value below this
#'   threshold is replaced by \code{NA}. Then set a maximum percent missing loci
#'   per sample: any sample with greater than this percentage missing loci is
#'   dropped. Note that threshold values can be explored without applying any
#'   filtering using the \code{explore_filter_coverage_samples()} function.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param min_coverage any coverage value \code{< min_coverage} is replaced by \code{NA}.
#' @param max_missing any sample with \code{> max_missing} values equal \code{NA} are dropped.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_coverage_samples <- function(x,
                                    min_coverage = 5,
                                    max_missing = 50,
                                    description = "filter samples based on coverage") {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_missing)
  assert_bounded(max_missing, right = 100)
  
  # replace low coverage with NA and drop some samples
  w <- which(x$coverage < min_coverage, arr.ind = TRUE)
  x$coverage[w] <- NA
  x$counts[w] <- NA
  percent_missing <- rowMeans(is.na(x$coverage)) * 100
  x$coverage <- x$coverage[percent_missing <= max_missing,]
  x$counts <- x$counts[percent_missing <= max_missing,]
  x$samples <- x$samples[percent_missing <= max_missing,]
  
  # record filter
  function_call <- paste0(deparse(match.call()), collapse = "")
  x$filter_history <- rbind(x$filter_history, list(description = description,
                                                   samples = nrow(x$coverage),
                                                   loci = ncol(x$coverage),
                                                   n_missing = sum(is.na(x$coverage)),
                                                   prop_missing = mean(is.na(x$coverage)),
                                                   function_call = function_call))
  
  return(x)
}

#------------------------------------------------
#' @title Explore locus coverage prior to filtering
#'
#' @description See \code{filter_coverage_loci()} function first. Explore values
#'   of \code{min_coverage} and \code{max_missing} without applying any
#'   filtering.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param min_coverage any coverage value \code{< min_coverage} is replaced by \code{NA}.
#' @param max_missing any locus with \code{> max_missing} values equal \code{NA} are dropped.
#' @param breaks number of breaks spanning the range \code{[0,100]}.
#'
#' @export

explore_filter_coverage_loci <- function(x,
                                         min_coverage = 5,
                                         max_missing = 50,
                                         breaks = 100) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_missing)
  assert_bounded(max_missing, right = 100)
  
  # get percent missing per locus
  x$coverage[x$coverage < min_coverage] <- NA
  percent_missing <- colMeans(is.na(x$coverage)) * 100
  percent_drop <- round(mean(percent_missing > max_missing)*100, digits = 2)
  
  # construct main title
  main_title <- paste0("min_coverage = ", min_coverage,
                       "\nmax_missing = ", max_missing,
                       "\nloci dropped = ", percent_drop, "%")
  
  # produce plot
  plot1 <- ggplot2::ggplot(data = data.frame(x = percent_missing)) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_histogram(ggplot2::aes(x = x), fill = "#4575B4", breaks = seq(0,100,l = breaks))
  plot1 <- plot1 + ggplot2::geom_vline(xintercept = max_missing, linetype = "dashed")
  plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(0,100))
  plot1 <- plot1 + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0))
  plot1 <- plot1 + ggplot2::xlab("% missing samples per locus")
  plot1 <- plot1 + ggplot2::ggtitle(main_title)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Filter loci based on coverage
#'
#' @description Set a coverage threshold: any coverage value below this
#'   threshold is replaced by \code{NA}. Then set a maximum percent missing samples
#'   per locus: any locus with greater than this percentage missing samples is
#'   dropped. Note that threshold values can be explored without applying any
#'   filtering using the \code{explore_filter_coverage_loci()} function.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param min_coverage any coverage value \code{< min_coverage} is replaced by \code{NA}.
#' @param max_missing any sample with \code{> max_missing} values equal \code{NA} are dropped.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_coverage_loci <- function(x,
                                 min_coverage = 5,
                                 max_missing = 50,
                                 description = "filter loci based on coverage") {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_missing)
  assert_bounded(max_missing, right = 100)
  
  # replace low coverage with NA and drop some samples
  w <- which(x$coverage < min_coverage, arr.ind = TRUE)
  x$coverage[w] <- NA
  x$counts[w] <- NA
  percent_missing <- colMeans(is.na(x$coverage)) * 100
  x$coverage <- x$coverage[,percent_missing <= max_missing]
  x$counts <- x$counts[,percent_missing <= max_missing]
  x$loci <- x$loci[percent_missing <= max_missing,]
  
  # record filter
  function_call <- paste0(deparse(match.call()), collapse = "")
  x$filter_history <- rbind(x$filter_history, list(description = description,
                                                   samples = nrow(x$coverage),
                                                   loci = ncol(x$coverage),
                                                   n_missing = sum(is.na(x$coverage)),
                                                   prop_missing = mean(is.na(x$coverage)),
                                                   function_call = function_call))
  
  return(x)
}

#------------------------------------------------
#' @title Plot coverage matrix
#'
#' @description Plot matrix of coverage over all samples and loci.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#'
#' @export

plot_coverage <- function(x) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  
  # ggplot raster
  plot1 <- ggplot2::ggplot(reshape2::melt(log(x$coverage)))
  plot1 <- plot1 + ggplot2::geom_raster(aes(x = Var1, y = Var2, fill = value))
  plot1 <- plot1 + ggplot2::scale_fill_viridis_c(option = "plasma", name = "log(coverage)")
  plot1 <- plot1 + ggplot2::theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  plot1 <- plot1 + ggplot2::xlab("locus") + ggplot2::ylab("sample")
  
  return(plot1)
}

#------------------------------------------------
#' @title Get within-sample allele frequencies
#'
#' @description Get within-sample allele frequencies from coverage and count
#'   data. Missing values can optionally be imputed by applying a summary
#'   function to the non NA values at each locus. The default summary function
#'   takes the mean of the non NA values.
#'
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param impute whether to impute missing values.
#' @param FUN function used to impute missing values. Default = `median`
#' @param ... other arguments to pass to \code{FUN}.
#' 
#' @export

get_wsaf <- function(x, impute = TRUE, FUN = median, ...) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_single_logical(impute)
  
  # get within-sample allele frequencies
  wsaf <- x$counts/x$coverage
  
  # impute missing values over loci
  if (impute) {
    locus_impute <- apply(wsaf, 2, FUN, na.rm = TRUE, ...)
    locus_impute <- outer(rep(1, nrow(wsaf)), locus_impute)
    wsaf[is.na(wsaf)] <- locus_impute[is.na(wsaf)]
  }
  
  return(wsaf)
}

#------------------------------------------------
#' @title Get genomic distance between samples
#'
#' @description Get genomic distance between samples using a distance metric
#'   that allows for mixed infections and takes account of linkage (see
#'   references for details).
#' 
#' @references MalariaGEN Plasmodium falciparum Community Project. "Genomic
#'   epidemiology of artemisinin resistant malaria". eLIFE (2016).
#' 
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param cutoff when calculating weights, correlations below this value are
#'   ignored (see references).
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console.
#'
#' @importFrom stats prcomp
#' @export

get_genomic_distance <- function(x, cutoff = 0.1, report_progress = TRUE) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_in(c("CHROM", "POS"), names(x$loci), message = "the `loci` component of the data must contain columns `CHROM` and `POS`")
  
  # get basic quantities
  nsamp <- nrow(x$samples)
  chrom_levels <- levels(x$loci$CHROM)
  nchrom <- length(chrom_levels)
  
  # get within-sample allele frequencies
  wsaf <- get_wsaf(x, impute = FALSE)
  
  # calculate weights
  weight <- NULL
  for (i in 1:nchrom) {
    
    # get correlation matrix between SNPs
    w <- which(x$loci$CHROM == chrom_levels[i])
    cormat <- suppressWarnings(cor(wsaf[,w], use = "pairwise.complete.obs"))
    cormat_na <- apply(cormat, 1, function(x) all(is.na(x)))
    cormat[cormat < cutoff] <- 0
    
    # calculate weights
    tmp <- rowSums(cormat^2, na.rm = TRUE)
    tmp[cormat_na] <- NA
    weight <- c(weight, 1/tmp)
  }
  
  # initialise progress bar
  if (report_progress) {
    pbar <- txtProgressBar(0, nsamp, style = 3)
  }
  
  # calculate distance matrix
  dist_mat <- matrix(NA, nsamp, nsamp)
  for (i in 1:(nsamp-1)) {
    
    # calculate distance between sample i and all others
    dab <- sweep(1-wsaf[-(1:i),,drop = FALSE], 2, wsaf[i,], "*") + sweep(wsaf[-(1:i),,drop = FALSE], 2, 1-wsaf[i,], "*")
    dab_weighted <- sweep(dab, 2, weight, "*")
    dist_mat[i,-(1:i)] <- rowSums(dab_weighted, na.rm = TRUE)
    
    # update progress bar
    if (report_progress) {
      setTxtProgressBar(pbar, i)
    }
  }
  
  # return distance matrix
  invisible(dist_mat)
}

#------------------------------------------------
#' @title PCA of within-sample allele frequencies
#'
#' @description Conduct principal components analysis (PCA) on a matrix of
#'   within-sample allele frequencies (WSAF). Missing values must have been
#'   already imputed. Output includes the raw components, the variance in the
#'   data explained by each component, and the loadings of each component also
#'   returned.
#'  
#' @param x a matrix of within-sample allele frequencies, as produced by the
#'   function \code{get_wsaf()}.
#'
#' @return Invisibly returns a list of class `prcomp` with the following 
#'   components
#'   \itemize{
#'       \item{"sdev"}{ the standard deviations of the principal components
#'       (i.e., the square roots of the eigenvalues of the
#'       covariance/correlation matrix, though the calculation is actually done
#'       with the singular values of the data matrix).} \item{"rotation"}{ the
#'       matrix of variable loadings (i.e., a matrix whose columns contain the
#'       eigenvectors). The function \code{princomp} returns this in the element
#'       \code{loadings}.} \item{"center, scale"}{ the centering and scaling
#'       used.} \item{"x"}{ the value of the rotated data (the centred data
#'       multiplied by the rotation matrix). Hence, \code{cov(x)} is the
#'       diagonal matrix \code{diag(sdev^2)}.}
#'       \item{"var"}{ the variance in the data explained by each component.}
#'       \item{"loadings"}{ the loadings of each component.}
#'       }
#' @importFrom stats prcomp
#' @export

pca_wsaf <- function(x) {
  
  # check inputs
  assert_matrix(x)
  assert_numeric(x)
  if (any(is.na(x))) {
    stop("input matrix cannot contain missing values. If produced using the function get_wsaf() then ensure that imputation is turned on")
  }
  
  # compute PCA
  pca <- prcomp(x)
  
  # compute variance explained
  pca$var <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2) * 100
  
  # compute loadings from rotations
  pca$loadings <- abs(pca$rotation)
  pca$loadings <- sweep(pca$loadings, 2, colSums(pca$loadings), "/")
  
  # return
  return(pca)
}

#------------------------------------------------
#' @title Plot variance explained by PCA components
#'
#' @description Plot the variance explained by each PCA component. The number of
#'   components shown is controlled by \code{num_components}, with up to the
#'   first 10 componenets shown by default. If less than the requested number of
#'   components exist, then all the components will be shown.
#'
#' @param pca output of \code{pca_wsaf()} function.
#' @param num_components maximum components to be shown.
#'
#' @importFrom plotly plot_ly layout
#' @export

plot_pca_variance <- function(pca, num_components = 10) {
  
  # potential components
  nc <- min(length(colnames(pca$x)), num_components)
  
  # x has to be factored in order to preserve PC order
  x <- factor(colnames(pca$x)[seq_len(nc)],
              levels = colnames(pca$x)[seq_len(nc)])
  
  # create plot object
  out <- plot_ly(x = x, y = pca$var[seq_len(nc)],
                 type = "bar", width = 500, height = 500) %>%
    plotly::layout(title = "Percentage of total variance explaned by PCA",
                   xaxis = list(title = "Principal Component"),
                   yaxis = list(title = "Percentage of variance explained"))
  
  # render and return invisibly
  print(out)
  invisible(out)
}

#------------------------------------------------
#' @title Plot PCA
#'
#' @description Plots either the first 2 or 3 principal components.
#'
#' @param pca output of \code{pca_wsaf()} function.
#' @param num_components numeric for number of components used. Default = 2.
#' @param col vector by which samples are coloured.
#' @param col_palette vector of colours for each group.
#'
#' @importFrom plotly plot_ly
#' @importFrom RColorBrewer brewer.pal
#' @export

plot_pca <- function(pca, num_components = 2, col = NULL, col_palette = NULL) {
  
  # check inputs
  assert_custom_class(pca, "prcomp")
  assert_in(num_components, c(2,3))
  if (is.null(col)) {
    col <- rep(1, nrow(pca$x))
  }
  assert_vector(col)
  assert_length(col, nrow(pca$x))
  if (is.null(col_palette)) {
    col_palette <- suppressWarnings(brewer.pal(3, "Set1"))
  }
  
  # check num_components
  nc <- min(ncol(pca$x), num_components)
  if (nc != num_components) {
    message("Using all the components (", nc, ") that are available.")
  }
  
  if (num_components == 2) {
    # scatterplot of first 2 principal components
    # 2D scatter
    plot1 <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2,
                     color = col, type = "scatter", colors = col_palette,
                     mode = "markers", marker = list(size = 5))
  }
  
  if (num_components == 3) {
    # scatterplot of first 3 principal components
    # 3D scatter
    plot1 <- plot_ly(as.data.frame(pca$x[, 1:3]), x = ~PC1, y = ~PC2, z = ~PC3,
                     color = col, type = "scatter3d", colors = col_palette,
                     mode = "markers", marker = list(size = 3))
  }
  
  # render and return invisibly
  print(plot1)
  invisible(plot1)
}

#------------------------------------------------
#' @title PCoA of genomic distances between samples
#'
#' @description Conduct principal coordinate analysis (PCoA) on a matrix of
#'   genomic distances.
#'  
#' @param x matrix of genomic distances, as produced by the function
#'   \code{get_genomic_distance()}.
#'
#' @importFrom ape pcoa
#' @export

pcoa_genomic_distance <- function(x) {
  
  # check inputs
  assert_matrix(x)
  assert_numeric(x)
  
  # mask the diagonal and reflect
  x[is.na(x)] <- 0
  x <- x + t(x)
  
  # compute PCoA
  ret <- ape::pcoa(x)
  
  # return PCoA
  invisible(ret)
}

#------------------------------------------------
#' @title Plot PCoA
#'
#' @description Plots either the first 2 or 3 vectors of PCoA.
#'
#' @param pcoa object of class "pcoa", as produced by \code{pcoa_wsaf()} function.
#' @param num_components numeric for number of components used. Default = 2.
#' @param col vector by which samples are coloured.
#' @param col_palette vector of colours for each group.
#'
#' @importFrom plotly plot_ly
#' @importFrom RColorBrewer brewer.pal
#' @export

plot_pcoa <- function(pcoa, num_components = 2, col = NULL, col_palette = NULL) {
  
  # check inputs
  assert_custom_class(pcoa, "pcoa")
  assert_in(num_components, c(2,3))
  if (is.null(col)) {
    col <- rep(1, nrow(pcoa$vectors))
  }
  assert_vector(col)
  assert_length(col, nrow(pcoa$vectors))
  if (is.null(col_palette)) {
    col_palette <- suppressWarnings(brewer.pal(3, "Set1"))
  }
  
  df <- data.frame(pcoa$vectors[,1:3])
  names(df) <- c("PC1", "PC2", "PC3")
  if (num_components == 2) {
    # scatterplot of first 2 principal components
    # 2D scatter
    plot1 <- plot_ly(df, x = ~PC1, y = ~PC2,
                     color = col, type = "scatter", colors = col_palette,
                     mode = "markers", marker = list(size = 5))
  }
  
  if (num_components == 3) {
    # scatterplot of first 3 principal components
    # 3D scatter
    plot1 <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
                     color = col, type = "scatter3d", colors = col_palette,
                     mode = "markers", marker = list(size = 3))
  }
  
  # render and return invisibly
  print(plot1)
  invisible(plot1)
}

#------------------------------------------------
#' @title Estimate pairwise inbreeding coefficient F by maximum likelihood
#'
#' @description Estimates the inbreeding coefficient between all pairs of
#'   samples by maximum likelihood.
#'
#' @details The probability of seeing the same or different alleles at a locus
#'   can be written in terms of the global allele frequency p and the inbreeding
#'   coefficient f, for example the probability of seeing the same REF allele is
#'   \code{(1-f)*p^2 + f*p}. This formula can be multiplied over all loci to
#'   arrive at the overall likelihood of each value of f, which can then be
#'   chosen by maximum likelihood. This function carries out this comparison
#'   between all pairwise samples, passed in as a matrix. The formula above only
#'   applies when comparing homozygous calls - for homo/het or het/het
#'   comparisons we can either ignore these loci (the default) or convert hets
#'   to homo by calling the major allele at every locus.
#'
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param f values of f that are explored.
#' @param ignore_het whether to ignore heterzygous comparisons, or alternatively
#'   call the major allele at every locus (see details).
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console.
#'
#' @export

inbreeding_mle <- function(x, f = seq(0,1,l=11), ignore_het = TRUE, report_progress = TRUE) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_vector(f)
  assert_bounded(f)
  assert_single_logical(ignore_het)
  assert_single_logical(report_progress)
  
  # get within-sample allele frequencies
  wsaf <- get_wsaf(x, impute = FALSE)
  
  # get global allele frequencies
  p <- colMeans(wsaf, na.rm = TRUE)
  
  # process hets
  if (ignore_het) {
    wsaf[wsaf != 0 & wsaf != 1] <- NA
  } else {
    wsaf <- round(wsaf)
  }
  
  # convert NA to -1 before passing to C++
  wsaf[is.na(wsaf)] <- -1
  
  # create progress bars
  pb <- txtProgressBar(min = 0, max = nrow(wsaf)-1, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # run efficient C++ function
  args <- list(x = mat_to_rcpp(wsaf), f = f, p = p, report_progress = report_progress)
  args_functions <- list(update_progress = update_progress)
  output_raw <- inbreeding_mle_cpp(args, args_functions, args_progress)
  
  # process output
  ret <- rcpp_to_mat(output_raw$ret)
  ret[row(ret) >= col(ret)] <- NA
  
  return(ret)
}

#------------------------------------------------
#' @title Simulate biallelic data
#'
#' @description Simulate biallelic data.
#'
#' @details TODO
#'
#' @param COI complexity of infection.
#' @param PLAF vector of population-level allele frequencies.
#' @param coverage coverage at each locus. If a single value then the same
#'   coverage is applied over all loci.
#' @param alpha shape parameter of the symmetric Dirichlet prior on strain
#'   proportions.
#' @param overdispersion the extent to which counts are over-dispersed relative
#'   to the binomial distribution. Counts are Beta-binomially distributed, with
#'   the beta distribution having shape parameters \code{p/overdispersion} and
#'   \code{(1-p)/overdispersion}.
#' @param epsilon the probability of a single read being mis-called as the other
#'   allele. Applies in both directions.
#'
#' @export

sim_biallelic <- function(COI = 3,
                          PLAF = seq(0,0.5,0.01),
                          coverage = 100,
                          alpha = 1,
                          overdispersion = 0,
                          epsilon = 0) {
  
  # check inputs
  assert_single_pos_int(COI)
  assert_vector(PLAF)
  assert_bounded(PLAF)
  L <- length(PLAF)
  if (length(coverage) == 1) {
    coverage <- rep(coverage, L)
  }
  assert_vector(coverage)
  assert_pos_int(coverage)
  assert_same_length(PLAF, coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion, zero_allowed = TRUE)
  assert_single_pos(epsilon, zero_allowed = TRUE)
  assert_bounded(epsilon)
  
  # generate strain proportions
  w <- rdirichlet(rep(alpha, COI))
  
  # generate true WSAF levels by summing binomial draws over strain proportions
  m <- mapply(function(x) rbinom(COI, 1, x), x = PLAF)
  p_levels <- colSums(sweep(m, 1, w, "*"))
  
  # add in genotyping error
  p_error <- p_levels*(1-epsilon) + (1-p_levels)*epsilon
  
  # draw read counts, taking into account overdispersion
  if (overdispersion == 0) {
    counts <- rbinom(L, size = coverage, prob = p_error)
  } else {
    counts <- rbetabinom(L, k = coverage, alpha = p_error/overdispersion, beta = (1-p_error)/overdispersion)
  }
  
  # return list
  ret <- list(COI = COI,
              strain_proportions = w,
              data = data.frame(PLAF = PLAF,
                                coverage = coverage,
                                counts = counts,
                                WSAF = counts/coverage))
  return(ret)
}

#------------------------------------------------
#' @title Dummy function
#'
#' @description Dummy function
#' 
#' @param x nubers to square.
#'
#' @export

dummy1 <- function(x = 1:5) {
  
  # print message to console
  message("running R dummy1 function")
  
  # get arguments in list form
  args <- list(x = x)
  
  # run C++ function with these arguments
  output_raw <- dummy1_cpp(args)
  
  # some optional processing of output
  message("processing output")
  ret <- output_raw$x_squared
  
  # return
  return(ret)
}


