
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
#' @title Ordinary print function for class mipanalyzer_multiallelic
#'
#' @description Calling \code{print()} on an object of class
#'   \code{mipanalyzer_multiallelic} results in custom output. This function therefore
#'   stands in for the base \code{print()} function, and is equivalent to
#'   running \code{print(unclass(x))}.
#'
#' @param x object of class \code{mipanalyzer_multiallelic}
#' @param ... other arguments passed to \code{print()}
#'
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_multiallelic")
  
  # print un-classed object
  print(unclass(x), ...)
  
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
