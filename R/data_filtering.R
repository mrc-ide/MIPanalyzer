#------------------------------------------------
#' @title Filter out some samples
#'
#' @description Filter out some samples.
#' 
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param sample_filter boolean vector specifying whether to keep (\code{TRUE})
#'   or drop (\code{FALSE}) each sample.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_samples <- function(x, sample_filter, description = "") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_vector(sample_filter)
  assert_logical(sample_filter)
  assert_eq(length(sample_filter), nrow(x$samples))
  
  # apply filter
  x$coverage <- x$coverage[sample_filter,,drop = FALSE]
  switch (class(x),
          "mipanalyzer_biallelic" = x$counts <- x$counts[sample_filter,,drop = FALSE],
          "mipanalyzer_multiallelic" = x$counts <- x$counts[,sample_filter,,drop = FALSE]
  )
  x$samples <- x$samples[sample_filter,,drop = FALSE]
  
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
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param locus_filter boolean vector specifying whether to keep (\code{TRUE})
#'   or drop (\code{FALSE}) each locus.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_loci <- function(x, locus_filter, description = "") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_vector(locus_filter)
  assert_logical(locus_filter)
  assert_eq(length(locus_filter), nrow(x$loci))
  
  # apply filter
  x$coverage <- x$coverage[,locus_filter,drop = FALSE]
  switch (class(x),
          "mipanalyzer_biallelic" = x$counts <- x$counts[,locus_filter,drop = FALSE],
          "mipanalyzer_multiallelic" = x$counts <- x$counts[,,locus_filter,drop = FALSE]
  )
  x$loci <- x$loci[locus_filter,,drop = FALSE]
  
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
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @export

filter_overcounts <- function(x, description = "replace overcounts with NA") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  
  # replace over-counts with NA
  
  # switch based on data type
  if (is(x, "mipanalyzer_biallelic")) {
    
    w <- which(x$counts > x$coverage, arr.ind = TRUE)
    x$coverage[w] <- NA
    x$counts[w] <- NA
    
  } else {
    
    for (i in seq_along(x$counts)[1]) {
      
      w <- which(x$counts[i,,] > x$coverage, arr.ind = TRUE)
      x$coverage[w] <- NA
      x$counts[i,,][w] <- NA
      
    }
  }
  
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
#' @title Filter alleles based on raw counts
#'
#' @description Drop any allele for which the number of read counts is below a
#'   given threshold. Coverage is adjusted to account for dropped reads.
#' 
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param count_min alleles with fewer than this many counts are dropped.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @importFrom methods is
#' @export

filter_counts <- function(x, count_min = 2, description = "filter individual allele counts") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_single_pos_int(count_min, zero_allowed = FALSE)
  
  # switch based on data type
  if (is(x, "mipanalyzer_biallelic")) {
    
    # drop alleles below threshold
    w <- which(!is.na(x$counts) & x$counts < count_min, arr.ind = TRUE)
    x$coverage[w] <- x$coverage[w] - x$counts[w]
    x$counts[w] <- 0
    w <- which(!is.na(x$counts) & (x$coverage - x$counts) < count_min, arr.ind = TRUE)
    x$coverage[w] <- x$counts[w]
    x$coverage[x$coverage == 0] <- NA
    
  } else {
    
    # drop alleles below threshold
    x$counts[!is.na(x$counts) & x$counts < count_min] <- 0
    x$coverage <- colSums(x$counts, na.rm = TRUE)
    x$coverage[x$coverage == 0] <- NA
  }
  
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
#' @title Filter alleles based on within-sample allele frequencies
#'
#' @description Drop any allele for which the within-sample allele frequency
#'   (WSAF) is below a givin threshold. Thresholds apply in both directions, for
#'   example if \code{wsaf_min = 0.01} then alleles with a WSAF less than 0.01
#'   *or* greater than 0.99 will be rounded to 0 or 1, respectively. Coverage is
#'   adjusted to account for dropped reads.
#' 
#' @param x object of class \code{mipanalyzer_multiallelic}.
#' @param wsaf_min alleles with counts that make a WSAF less than this threshold
#'   are dropped.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @importFrom methods is
#' @export

filter_wsaf <- function(x, wsaf_min = 0.01, description = "filter individual allele WSAF") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_single_bounded(wsaf_min)
  
  # get WSAF
  wsaf <- get_wsaf(x, impute = FALSE)
  
  # switch based on data type
  if (is(x, "mipanalyzer_biallelic")) {
    
    # drop alleles below threshold
    w <- which(!is.na(wsaf) & wsaf < wsaf_min, arr.ind = TRUE)
    x$coverage[w] <- x$coverage[w] - x$counts[w]
    x$counts[w] <- 0
    w <- which(!is.na(wsaf) & (1.0 - wsaf) < wsaf_min, arr.ind = TRUE)
    x$coverage[w] <- x$counts[w]
    x$coverage[x$coverage == 0] <- NA
    x$counts[w] <- NA
    
  } else {
    
    # drop alleles below threshold
    for (i in 1:4) {
      x$counts[i,,][!is.na(wsaf[i,,]) & wsaf[i,,] < wsaf_min] <- NA
    }
    x$coverage <- colSums(x$counts, na.rm = TRUE)
    x$coverage[x$coverage == 0] <- NA
  }
  
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
#' @title Filter loci to drop invariant sites
#'
#' @description Filter loci to drop invariant sites.
#' 
#' @param x object of class \code{mipanalyzer_multiallelic}.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @importFrom methods is
#' @export

filter_loci_invariant <- function(x, description = "filter loci to drop invariant sites") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  
  # get WSAF
  wsaf <- get_wsaf(x, impute = FALSE)
  
  # switch based on data type
  if (is(x, "mipanalyzer_biallelic")) {
    
    # identify invariant sites
    invariant <- apply(wsaf, 2, function(y) {
      all(y[!is.na(y)] == 1) | all(y[!is.na(y)] == 0)
    })
    
  } else {
    
    # identify invariant sites
    invariant <- 1
    for (i in 1:dim(wsaf)[1]) {
      invariant_i <- apply(wsaf[i,,], 2, function(y) {
        all(y[!is.na(y)] == 1) | all(y[!is.na(y)] == 0)
      })
      invariant <- invariant*invariant_i
    }
  }
  
  # drop invariant sites
  ret <- filter_loci(x, !invariant, description = description)
  
  return(ret)
}

#------------------------------------------------
#' @title Explore sample coverage prior to filtering
#'
#' @description Explore what effect the \code{filter_coverage_samples()}
#'   function will have on the data without actually applying any filters. Can
#'   be used to set coverage thresholds.
#' 
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param min_coverage the coverage threshold below which data is deemed to be
#'   low-coverage.
#' @param max_low_coverage (percentage). Samples are not allowed to contain more
#'   than this many low-coverage loci. In the \code{filter_coverage_samples()}
#'   function, any sample with more than \code{max_low_coverage} low-coverage
#'   loci will be dropped.
#' @param breaks number of breaks spanning the range \code{[0,100]}.
#'
#' @export

explore_filter_coverage_samples <- function(x,
                                            min_coverage = 5,
                                            max_low_coverage = 50,
                                            breaks = 100) {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_low_coverage)
  assert_bounded(max_low_coverage, right = 100)
  
  # get percent low-coverage loci per sample
  percent_low_coverage <- rowMeans(x$coverage < min_coverage | is.na(x$coverage)) * 100
  percent_drop <- round(mean(percent_low_coverage > max_low_coverage)*100, digits = 2)
  
  # construct main title
  main_title <- paste0("min_coverage = ", min_coverage,
                       "\nmax_low_coverage = ", max_low_coverage, "%",
                       "\nsamples dropped = ", percent_drop, "%")
  
  # produce plot
  plot1 <- ggplot2::ggplot(data = data.frame(x = percent_low_coverage)) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_histogram(ggplot2::aes(x = x), fill = "#4575B4", breaks = seq(0,100,l = breaks))
  plot1 <- plot1 + ggplot2::geom_vline(xintercept = max_low_coverage, linetype = "dashed")
  plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(0,100))
  plot1 <- plot1 + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0))
  plot1 <- plot1 + ggplot2::xlab("% low-coverage loci per sample")
  plot1 <- plot1 + ggplot2::ggtitle(main_title)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Filter samples based on coverage
#'
#' @description Set a coverage threshold: any coverage value below this
#'   threshold is deemed to be low-coverage. Then set a maximum percent
#'   low-coverage loci per sample: any sample with greater than this percentage
#'   low-coverage loci is dropped. Note that threshold values can be explored
#'   without applying any filtering using the
#'   \code{explore_filter_coverage_samples()} function.
#' 
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param min_coverage the coverage threshold below which data is deemed to be
#'   low-coverage.
#' @param max_low_coverage any sample with more than \code{max_low_coverage}
#'   percent of low-coverage loci will be dropped.
#' @param replace_low_coverage (Boolean). If \code{TRUE} then any remaining
#'   low-coverage loci will be replaced with \code{NA}.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @importFrom methods is
#' @export

filter_coverage_samples <- function(x,
                                    min_coverage = 5,
                                    max_low_coverage = 50,
                                    replace_low_coverage = FALSE,
                                    description = "filter samples based on coverage") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_low_coverage)
  assert_bounded(max_low_coverage, right = 100)
  assert_single_logical(replace_low_coverage)
  
  # get percent low-coverage loci per sample
  percent_low_coverage <- rowMeans(x$coverage < min_coverage | is.na(x$coverage)) * 100
  
  # drop samples with too many low-coverage loci
  x$coverage <- x$coverage[percent_low_coverage <= max_low_coverage,,drop = FALSE]
  switch (class(x),
          "mipanalyzer_biallelic" = x$counts <- x$counts[percent_low_coverage <= max_low_coverage,,drop = FALSE],
          "mipanalyzer_multiallelic" = x$counts <- x$counts[,percent_low_coverage <= max_low_coverage,,drop = FALSE]
  )
  x$samples <- x$samples[percent_low_coverage <= max_low_coverage,,drop = FALSE]
  
  # replace low-coverage with NA
  if (replace_low_coverage) {
    w <- which(x$coverage < min_coverage, arr.ind = TRUE)
    x$coverage[w] <- NA
    if (is(x, "mipanalyzer_multiallelic")) {
      for (i in 1:4) {
        x$counts[i,,][w] <- NA
      }
    } else {
      x$counts[w] <- NA
    }
  }
  
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
#' @description Explore what effect the \code{filter_coverage_loci()}
#'   function will have on the data without actually applying any filters. Can
#'   be used to set coverage thresholds.
#' 
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param min_coverage the coverage threshold below which data is deemed to be
#'   low-coverage.
#' @param max_low_coverage (percentage). Loci are not allowed to contain more
#'   than this many low-coverage samples. In the \code{filter_coverage_loci()}
#'   function, any locus with more than \code{max_low_coverage} low-coverage
#'   samples will be dropped.
#' @param breaks number of breaks spanning the range \code{[0,100]}.
#'
#' @export

explore_filter_coverage_loci <- function(x,
                                         min_coverage = 5,
                                         max_low_coverage = 50,
                                         breaks = 100) {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_low_coverage)
  assert_bounded(max_low_coverage, right = 100)
  
  
  # get percent low-coverage loci per sample
  percent_low_coverage <- colMeans(x$coverage < min_coverage | is.na(x$coverage)) * 100
  percent_drop <- round(mean(percent_low_coverage > max_low_coverage)*100, digits = 2)
  
  # construct main title
  main_title <- paste0("min_coverage = ", min_coverage,
                       "\nmax_low_coverage = ", max_low_coverage, "%",
                       "\nloci dropped = ", percent_drop, "%")
  
  # produce plot
  plot1 <- ggplot2::ggplot(data = data.frame(x = percent_low_coverage)) + ggplot2::theme_bw()
  plot1 <- plot1 + ggplot2::geom_histogram(ggplot2::aes(x = x), fill = "#4575B4", breaks = seq(0,100,l = breaks))
  plot1 <- plot1 + ggplot2::geom_vline(xintercept = max_low_coverage, linetype = "dashed")
  plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(0,100))
  plot1 <- plot1 + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0))
  plot1 <- plot1 + ggplot2::xlab("% low-coverage samples per locus")
  plot1 <- plot1 + ggplot2::ggtitle(main_title)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Filter loci based on coverage
#'
#' @description Set a coverage threshold: any coverage value below this
#'   threshold is deemed to be low-coverage. Then set a maximum percent
#'   low-coverage samples per locus: any locus with greater than this percentage
#'   low-coverage samples is dropped. Note that threshold values can be explored
#'   without applying any filtering using the
#'   \code{explore_filter_coverage_loci()} function.
#'
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#' @param min_coverage the coverage threshold below which data is deemed to be
#'   low-coverage.
#' @param max_low_coverage any locus with more than \code{max_low_coverage}
#'   percent of low-coverage samples will be dropped.
#' @param replace_low_coverage (Boolean). If \code{TRUE} then any remaining
#'   low-coverage loci will be replaced with \code{NA}.
#' @param description brief description of the filter, to be saved in the filter
#'   history.
#'
#' @importFrom methods is
#' @export

filter_coverage_loci <- function(x,
                                 min_coverage = 5,
                                 max_low_coverage = 50,
                                 replace_low_coverage = FALSE,
                                 description = "filter loci based on coverage") {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_single_pos_int(min_coverage)
  assert_single_pos(max_low_coverage)
  assert_bounded(max_low_coverage, right = 100)
  assert_single_logical(replace_low_coverage)
  
  # get percent low-coverage loci per sample
  percent_low_coverage <- colMeans(x$coverage < min_coverage | is.na(x$coverage)) * 100
  
  # drop loci with too many low-coverage samples
  x$coverage <- x$coverage[,percent_low_coverage <= max_low_coverage]
  switch (class(x),
          "mipanalyzer_biallelic" = x$counts <- x$counts[,percent_low_coverage <= max_low_coverage],
          "mipanalyzer_multiallelic" = x$counts <- x$counts[,,percent_low_coverage <= max_low_coverage,drop = FALSE]
  )
  x$loci <- x$loci[percent_low_coverage <= max_low_coverage,]
  
  # replace low-coverage with NA
  if (replace_low_coverage) {
    w <- which(x$coverage < min_coverage, arr.ind = TRUE)
    x$coverage[w] <- NA
    if (is(x, "mipanalyzer_multiallelic")) {
      for (i in 1:4) {
        x$counts[i,,][w] <- NA
      }
    } else {
      x$counts[w] <- NA
    }
  }
  
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
