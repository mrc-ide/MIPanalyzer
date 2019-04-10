
#------------------------------------------------
#' @title Custom print function for class mipanalyzer_multiallelic
#'   
#' @description Custom print function for class \code{mipanalyzer_multiallelic},
#'   printing a summary of the key elements (also equivalent to
#'   \code{summary(x)}). To do an ordinary \code{print()}, use the
#'   \code{print_full()} function.
#' 
#' @param x object of class \code{mipanalyzer_multiallelic}
#' @param ... other arguments (ignored)
#' 
#' @export

print.mipanalyzer_multiallelic <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class mipanalyzer_multiallelic
#'   
#' @description Custom summary function for class \code{mipanalyzer_multiallelic}.
#'   
#' @param object object of class \code{mipanalyzer_multiallelic}
#' @param ... other arguments (ignored)
#'   
#' @export

summary.mipanalyzer_multiallelic <- function(object, ...) {
  
  # summarise raw data
  n_samples <- nrow(object$coverage)
  n_loci <- ncol(object$coverage)
  
  message(sprintf("samples = %s\nloci = %s", n_samples, n_loci))
}

#------------------------------------------------
#' @title Determine if object is of class mipanalyzer_multiallelic
#'
#' @description Determine if object is of class \code{mipanalyzer_multiallelic}.
#'
#' @param x object of class \code{mipanalyzer_multiallelic}
#'
#' @export

is.mipanalyzer_multiallelic <- function(x) {
  inherits(x, "mipanalyzer_multiallelic")
}
