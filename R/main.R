
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
  counts <- masplit(counts_raw, record = 1)
  
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
#' @param FUN function used to impute missing values. Default = `mean`
#' @param ... other arguments to pass to \code{FUN}.
#' 
#' @export

get_wsaf <- function(x, impute = TRUE, FUN = mean, ...) {
  
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
#' @title PCA of within-sample allele frequencies
#'
#' @description Conduct principal components analysis (PCA) on input matrix of
#'   within-sample allele frequencies (WSAF), as returned by \code{get_WSAF()}
#'   function. Requires that any missing values have already been imputed.
#'   Ouptut includes the raw components, the variance in the data explained by
#'   each component, and the loadings of each component also returned.
#'  
#' @param x within-sample allele frequency matrix.
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
    stop("input matrix cannot contain any missing values")
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
#'
#' @importFrom plotly plot_ly
#' @importFrom RColorBrewer brewer.pal
#' @export

plot_pca <- function(pca, num_components = 2) {
  
  # check inputs
  assert_custom_class(pca, "prcomp")
  assert_in(num_components, c(2,3))
  
  # check num_components
  nc <- min(ncol(pca$x), num_components)
  if (nc != num_components) {
    message("Using all the components (", nc, ") that are available.")
  }
  
  # create color vector
  col_vec <- suppressWarnings(brewer.pal(3, "Set1"))
  
  if (nc == 2) {
    # scatterplot of first 2 principal components
    # 2D scatter
    plot1 <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2,
                     color = 1, type = "scatter", colors = col_vec,
                     mode = "markers", marker = list(size = 5))
  }
  
  if (nc == 3) {
    # scatterplot of first 3 principal components
    # 3D scatter
    plot1 <- plot_ly(as.data.frame(pca$x[, 1:3]), x = ~PC1, y = ~PC2, z = ~PC3,
                     color = 1, type = "scatter3d", colors = col_vec,
                     mode = "markers", marker = list(size = 3))
  }
  
  # render and return invisibly
  print(plot1)
  invisible(plot1)
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


