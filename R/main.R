
#------------------------------------------------
#' @useDynLib MIPanalyzer, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
#' @title Create mipanalyzer_data object
#'
#' @description Create mipanalyzer_data object.
#'
#' @param n number of samples.
#' @param L number of loci
#'
#' @export

mipanalyzer_data <- function(n = 50, L = 100) {
  
  # create dummy data
  ret <- list(coverage = array(0, dim = c(3, n, L)),
              counts = array(0, dim = c(3, n, L)),
              samples = data.frame(ID = 1:n),
              loci = data.frame(POS = 1:L))
  
  # return in mipanalyzer_data class
  class(ret) <- "mipanalyzer_data"
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


