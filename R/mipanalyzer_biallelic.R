
#------------------------------------------------
#' @title Custom print function for class mipanalyzer_biallelic
#'   
#' @description Custom print function for class \code{mipanalyzer_biallelic},
#'   printing a summary of the key elements (also equivalent to
#'   \code{summary(x)}). To do an ordinary \code{print()}, use the
#'   \code{print_full()} function.
#' 
#' @param x object of class \code{mipanalyzer_biallelic}
#' @param ... other arguments (ignored)
#' 
#' @export

print.mipanalyzer_biallelic <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class mipanalyzer_biallelic
#'   
#' @description Custom summary function for class \code{mipanalyzer_biallelic}.
#'   
#' @param object object of class \code{mipanalyzer_biallelic}
#' @param ... other arguments (ignored)
#'   
#' @export

summary.mipanalyzer_biallelic <- function(object, ...) {
  
  # summarise raw data
  n_samples <- nrow(object$coverage)
  n_loci <- ncol(object$coverage)
  
  message(sprintf("samples = %s\nloci = %s", n_samples, n_loci))
}

#------------------------------------------------
#' @title Determine if object is of class mipanalyzer_biallelic
#'
#' @description Determine if object is of class \code{mipanalyzer_biallelic}.
#'
#' @param x object of class \code{mipanalyzer_biallelic}
#'
#' @export

is.mipanalyzer_biallelic <- function(x) {
  inherits(x, "mipanalyzer_biallelic")
}
