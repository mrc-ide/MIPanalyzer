
#------------------------------------------------
#' @title Check that MIPanalyzer package has loaded successfully
#'
#' @description Simple function to check that MIPanalyzer package has loaded 
#'   successfully.
#'
#' @export

check_MIPanalyzer_loaded <- function() {
  message("MIPanalyzer  loaded successfully!")
}

#------------------------------------------------
#' @title Plot coverage matrix
#'
#' @description Plot matrix of coverage over all samples and loci.
#' 
#' @param x object of class \code{mipanalyzer_biallelic} or
#'   \code{mipanalyzer_multiallelic}.
#'
#' @import ggplot2
#' @export

plot_coverage <- function(x) {
  
  # avoid "no visible binding" notes
  Var1 <- Var2 <- value <- NULL
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  
  # ggplot raster
  plot1 <- ggplot(reshape2::melt(log(x$coverage)))
  plot1 <- plot1 + geom_raster(aes(x = Var1, y = Var2, fill = value))
  plot1 <- plot1 + scale_fill_viridis_c(option = "plasma", name = "log(coverage)")
  plot1 <- plot1 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  plot1 <- plot1 + xlab("sample") + ylab("locus")
  
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
#' @param FUN function used to impute missing values. Default = `median`
#' @param ... other arguments to pass to \code{FUN}.
#' 
#' @importFrom methods is
#' @export

get_wsaf <- function(x, impute = TRUE, FUN = median, ...) {
  
  # check inputs
  assert_custom_class(x, c("mipanalyzer_biallelic", "mipanalyzer_multiallelic"))
  assert_single_logical(impute)
  
  # switch based on class
  if (is(x, "mipanalyzer_biallelic")) {
    
    # get within-sample allele frequencies
    wsaf <- x$counts/x$coverage
    
    # impute missing values over loci
    if (impute) {
      locus_impute <- apply(wsaf, 2, FUN, na.rm = TRUE, ...)
      locus_impute <- outer(rep(1, nrow(wsaf)), locus_impute)
      wsaf[is.na(wsaf)] <- locus_impute[is.na(wsaf)]
    }
    
  } else {
    
    # get within-sample allele frequencies
    wsaf <- array(NA, dim = dim(x$counts))
    for (i in 1:4) {
      wsaf[i,,] <- x$counts[i,,]/x$coverage
    }
    
    # impute missing values over loci
    if (impute) {
      for (i in 1:4) {
        locus_impute <- apply(wsaf[i,,], 2, FUN, na.rm = TRUE, ...)
        locus_impute <- outer(rep(1, nrow(wsaf[i,,])), locus_impute)
        wsaf[i,,][is.na(wsaf[i,,])] <- locus_impute[is.na(wsaf[i,,])]
      }
    }
    
  }
  
  return(wsaf)
}

#------------------------------------------------
#' @title Produce ggplot map
#'
#' @description Produce ggplot map.
#'
#' @param x_limits longitude limits of map.
#' @param y_limits latitude limits of map.
#' @param col_country fill colour of countries.
#' @param col_country_border colour of country borders.
#' @param size_country_border size of country borders.
#' @param col_sea fill colour of sea.
#' @param resolution the resolution of the underlying map. Must be one of
#'   "coarse", "low", "less", "islands", "li", "high".
#'
#' @importFrom rworldmap getMap
#' @export

plot_map <- function(x_limits = c(12, 35),
                     y_limits = c(-13,5),
                     col_country = grey(0.3),
                     col_country_border = grey(0.5),
                     size_country_border = 0.5,
                     col_sea = grey(0.1),
                     resolution = "coarse") {
  
  # avoid "no visible binding" notes
  long <- lat <- group <- NULL
  
  # check inputs
  assert_vector(x_limits)
  assert_length(x_limits, 2)
  assert_numeric(x_limits)
  assert_vector(y_limits)
  assert_length(y_limits, 2)
  assert_numeric(y_limits)
  assert_in(resolution, c("coarse", "low", "less", "islands", "li", "high"))
  
  # load country shapefiles
  world_map <- getMap(resolution = resolution)
  
  # basic plot
  plot1 <- ggplot() + theme_bw() + theme(panel.background = element_rect(fill = col_sea),
                                         panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank())
  
  # add country borders
  plot1 <- plot1 + geom_polygon(aes(long, lat, group = group),
                                size = size_country_border, color = col_country_border,
                                fill = col_country, data = world_map)
  
  # limits and labels
  plot1 <- plot1 + xlab("longitude") + ylab("latitude")
  plot1 <- plot1 + coord_cartesian(xlim = x_limits, ylim = y_limits) 
  
  # return
  return(plot1)
}

#------------------------------------------------
#' @title Get dataframe of P.falciparum chromosome lengths
#'
#' @description Get dataframe of P.falciparum chromosome lengths
#'
#' @export

Pf_chrom_lengths <- function() {
  ret <- data.frame(chrom = 1:14,
                    length = c(643292, 947102, 1060087,
                               1204112, 1343552, 1418244,
                               1501717, 1419563, 1541723,
                               1687655, 2038337, 2271478,
                               2895605, 3291871))
  return(ret)
}

#------------------------------------------------
#' @title Ordinary print function for unclassed object
#'
#' @description Calling \code{print()} on an object of custom class, e.g.
#'   \code{mipanalyzer_biallelic}, results in custom output. This function
#'   therefore stands in for the base \code{print()} function, and is equivalent
#'   to running \code{print(unclass(x))}.
#'
#' @param x object of custom class (e.g. \code{mipanalyzer_biallelic}).
#' @param ... other arguments passed to \code{print()}.
#'
#' @export

print_full <- function(x, ...) {
  
  # print un-classed object
  print(unclass(x), ...)
  
  # return invisibly
  invisible(x)
}
