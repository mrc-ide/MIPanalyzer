#------------------------------------------------
#' @title Get great circle distance between spatial points
#'
#' @description Get great circle distance between spatial points.
#' 
#' @param lat vector of latitudes.
#' @param long vector of longitudes.
#'
#' @export

get_spatial_distance <- function(lat, long) {
  
  # check inputs
  assert_vector(lat)
  assert_numeric(lat)
  assert_vector(long)
  assert_numeric(long)
  
  # calculate distance matrix
  ret <- apply(cbind(lat, long), 1, function(y) {lonlat_to_bearing(lat, long, y[1], y[2])$gc_dist})
  diag(ret) <- 0
  
  return(ret)
}

#------------------------------------------------
#' @title Calculate great circle distance and bearing between coordinates
#'
#' @description Calculate great circle distance and bearing between spatial
#'   coordinates.
#'
#' @param origin_lon The origin longitude
#' @param origin_lat The origin latitude
#' @param dest_lon The destination longitude
#' @param dest_lat The destination latitude
#'
#' @export
#' @examples
#' # one degree longitude should equal approximately 111km at the equator
#' lonlat_to_bearing(0, 0, 1, 0)

lonlat_to_bearing <- function(origin_lon, origin_lat, dest_lon, dest_lat) {
  
  # check inputs
  assert_vector(origin_lon)
  assert_numeric(origin_lon)
  assert_vector(origin_lat)
  assert_numeric(origin_lat)
  assert_vector(dest_lon)
  assert_numeric(dest_lon)
  assert_vector(dest_lat)
  assert_numeric(dest_lat)
  
  # convert input arguments to radians
  origin_lon <- origin_lon*2*pi/360
  origin_lat <- origin_lat*2*pi/360
  dest_lon <- dest_lon*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  
  # get change in lon
  delta_lon <- dest_lon - origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or 
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  
  # convert bearing from radians to degrees measured clockwise from due north,
  # and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earth_rad <- 6371
  gc_dist <- earth_rad*gc_angle
  
  # return list
  ret <-list(bearing = bearing,
             gc_dist = gc_dist)
  return(ret)
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
  chrom_levels <- unique(x$loci$CHROM)
  nchrom <- length(chrom_levels)
  
  # get within-sample allele frequencies
  wsaf <- get_wsaf(x, impute = FALSE)
  
  # calculate weights
  weight <- NULL
  for (i in 1:nchrom) {
    
    # get correlation matrix between SNPs
    w <- which(x$loci$CHROM == chrom_levels[i])
    cormat <- suppressWarnings(cor(wsaf[,w,drop = FALSE], use = "pairwise.complete.obs"))
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
    dist_mat[i,-(1:i)] <- rowSums(dab_weighted, na.rm = TRUE) / apply(dab_weighted, 1, function(x) sum(!is.na(x)))
    
    # update progress bar
    if (report_progress) {
      setTxtProgressBar(pbar, i)
    }
  }
  
  # return distance matrix
  invisible(dist_mat)
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
#'   \eqn{(1-f)*p^2 + f*p}. This formula can be multiplied over all loci to
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

inbreeding_mle <- function(x, f = seq(0,1,l=11), ignore_het = FALSE, report_progress = TRUE) {
  
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
  ret_ml <- rcpp_to_mat(output_raw$ret_ml)
  ret_ml[row(ret_ml) >= col(ret_ml)] <- NA
  ret_all <- rcpp_to_array(output_raw$ret_all)
  
  # return list
  ret <- list(mle = ret_ml,
              loglike = ret_all)
  return(ret)
}


#------------------------------------------------
#' @title Get identity by mixture
#'
#' @description Get identity by "mixture distance". The mixture distance between
#'   two samples is the proportion of loci that have identical within-sample
#'   allele frequencies (WSAFs), or alternatively have WSAFs within a given
#'   tolerance. This extends the idea of identity by state (IBS) to continuous
#'   WSAFs rather than categorical genotypes.
#'
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param tol tolerance on mixture comparisons. Default = 0.
#' @param diagonal Should the diagonal of the distance matrix be changed to a
#'   given value. Default = NULL, which cause no changes.
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console.
#'
#' @export

get_IB_mixture <- function(x, tol = 0, diagonal = NULL, report_progress = TRUE) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_single_numeric(tol)
  if (!is.null(diagonal)) {
    assert_single_numeric(diagonal)
  }
  assert_single_logical(report_progress)
  
  # get basic quantities
  nsamp <- nrow(x$samples)
  
  # get within-sample allele frequencies
  wsaf <- get_wsaf(x, impute = FALSE)
  
  # initialise progress bar
  if (report_progress) {
    pbar <- txtProgressBar(0, nsamp, style = 3)
  }
  
  # compute all pairwise IBS
  ret <- matrix(NA, nsamp, nsamp)
  for (i in 1:nsamp) {
    if (tol == 0) {
      ret[i,] <- rowMeans(outer(rep(1,nsamp), wsaf[i,]) == wsaf, na.rm = TRUE)
    } else {
      ret[i,] <- rowMeans(abs(outer(rep(1,nsamp), wsaf[i,]) - wsaf) <= tol, na.rm = TRUE)  
    }
    # update progress bar
    if (report_progress) {
      setTxtProgressBar(pbar, i)
    }
  }
  
  if (!is.null(diagonal)) {
    ret[diag(ret)] <- diagonal
  }
  ret[lower.tri(ret)] <- NA
  
  # return distance matrix
  return(ret)
}


#------------------------------------------------
#' @title Get identity by state (IBS) distance
#'
#' @description Get identity by state (IBS) distance, computed as the proportion
#'   of sites that are identical between samples. If \code{ignore_het = TRUE}
#'   then heterozygous sites are ignored, otherwise the major strain is called
#'   at every locus.
#'
#' @param x object of class \code{mipanalyzer_biallelic}.
#' @param ignore_het whether to ignore heterozygous comparisons, or alternatively
#'   call the major allele at every locus (see details).
#' @param diagonal Should the diagonal of the distance matrix be changed to a
#'   given value. Default = \code{NULL}, which cause no changes.
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console.
#'
#' @export

get_IBS_distance <- function(x, ignore_het = TRUE, diagonal = NULL, report_progress = TRUE) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  assert_single_logical(ignore_het)
  assert_single_logical(report_progress)
  
  # get basic quantities
  nsamp <- nrow(x$samples)
  
  # get within-sample allele frequencies
  wsaf <- get_wsaf(x, impute = FALSE)
  
  # deal with heterozygous sites
  if (ignore_het) {
    wsaf[wsaf != 0 & wsaf != 1] <- NA
  } else {
    wsaf <- round(wsaf)
  }
  
  # initialise progress bar
  if (report_progress) {
    pbar <- txtProgressBar(0, nsamp, style = 3)
  }
  
  # compute all pairwise IBS
  ret <- matrix(NA, nsamp, nsamp)
  for (i in 1:nsamp) {
    ret[i,] <- rowMeans(outer(rep(1,nsamp), wsaf[i,]) == wsaf, na.rm = TRUE)
    
    # update progress bar
    if (report_progress) {
      setTxtProgressBar(pbar, i)
    }
  }
  
  if (!is.null(diagonal)) {
    ret[diag(ret)] <- diagonal
  }
  ret[lower.tri(ret)] <- NA
  
  
  # return distance matrix
  return(ret)
}

#------------------------------------------------
#' @title Plot a distance matrix
#'
#' @description Simple image plot of a matrix of pairwise distances.
#'
#' @param m square matrix of pairwise distances.
#' @param col_pal which viridis colour pallet to use. Options are "viridis",
#'   "plasma", "magma" or "inferno".
#'
#' @export

plot_distance <- function(m, col_pal = "plasma") {
  
  # avoid "no visible binding" notes
  x <- y <- value <- NULL
  
  # check inputs
  assert_square_matrix(m)
  assert_single_string(col_pal)
  assert_in(col_pal, c("viridis", "plasma", "magma", "inferno"))
  
  # get matrix into long-form data.frame
  n <- nrow(m)
  df_plot <- data.frame(x = 1:n,
                        y = rep(1:n, each = n),
                        value = as.vector(m))
  
  # produce plot
  df_plot |>
    ggplot() + theme_void() +
    geom_tile(aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_c(option = col_pal, na.value = "white")
  
}
