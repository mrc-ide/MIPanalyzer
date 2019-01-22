
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
#' @title Ordinary print function for class mipanalyzer_biallelic
#'
#' @description Calling \code{print()} on an object of class
#'   \code{mipanalyzer_biallelic} results in custom output. This function therefore
#'   stands in for the base \code{print()} function, and is equivalent to
#'   running \code{print(unclass(x))}.
#'
#' @param x object of class \code{mipanalyzer_biallelic}
#' @param ... other arguments passed to \code{print()}
#'
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_biallelic")
  
  # print un-classed object
  print(unclass(x), ...)
  
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
