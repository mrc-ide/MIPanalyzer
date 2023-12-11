#------------------------------------------------
#' @title MIPanalyzer
#'
#' @description This package can be used to read in raw molecular inversion
#'   probe (MIP) data from vcf into a format that is convenient to work with.
#'   Data can be filtered based on counts, frequencies, missingness or other
#'   criteria. Filtered data can be analysed by common methods including PCA and
#'   various pairwise genetic metrics, and can be visualised in multiple ways.
#'   This package is intended to evolve as new MIP analyses are needed, thereby
#'   making it easy to repeat common analyses as new data becomes available.
#'
#' @docType package
#' @name MIPanalyzer
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib MIPanalyzer, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("MIPanalyzer", libpath)
}
