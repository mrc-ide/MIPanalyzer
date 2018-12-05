
#------------------------------------------------
#' @title Custom print function for class mipanalyzer_data
#'   
#' @description Custom print function for class \code{mipanalyzer_data},
#'   printing a summary of the key elements (also equivalent to
#'   \code{summary(x)}). To do an ordinary \code{print()}, use the
#'   \code{print_full()} function.
#' 
#' @param x object of class \code{mipanalyzer_data}
#' @param ... other arguments (ignored)
#' 
#' @export

print.mipanalyzer_data <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Ordinary print function for class mipanalyzer_data
#'
#' @description Calling \code{print()} on an object of class
#'   \code{mipanalyzer_data} results in custom output. This function therefore
#'   stands in for the base \code{print()} function, and is equivalent to
#'   running \code{print(unclass(x))}.
#'
#' @param x object of class \code{mipanalyzer_data}
#' @param ... other arguments passed to \code{print()}
#'
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_custom_class(x, "mipanalyzer_data")
  
  # print un-classed object
  print(unclass(x), ...)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class mipanalyzer_data
#'   
#' @description Custom summary function for class \code{mipanalyzer_data}.
#'   
#' @param object object of class \code{mipanalyzer_data}
#' @param ... other arguments (ignored)
#'   
#' @export

summary.mipanalyzer_data <- function(object, ...) {
  
  # summarise raw data
  n_samples <- dim(object$coverage)[2]
  n_loci <- dim(object$coverage)[3]
  
  message(sprintf("samples = %s", n_samples))
  message(sprintf("loci = %s", n_loci))
}

#------------------------------------------------
#' @title Determine if object is of class mipanalyzer_data
#'
#' @description Determine if object is of class \code{mipanalyzer_data}.
#'
#' @param x object of class \code{mipanalyzer_data}
#'
#' @export

is.mipanalyzer_data <- function(x) {
  inherits(x, "mipanalyzer_data")
}
