#------------------------------------------------
#' @title PCA of within-sample allele frequencies
#' 
#' @description Conduct principal components analysis (PCA) on a matrix of
#'   within-sample allele frequencies (WSAF). Missing values must have been
#'   already imputed. Output includes the raw components, the variance in the
#'   data explained by each component, and the loadings of each component also
#'   returned.
#' 
#' @details Contributions of each variable are computed from the loading values
#'   (stored as "rotation" within the \code{prcomp} object). The percent
#'   contribution of a variable is defined as the absolute loading value for
#'   this variable, divided by the sum of loadings over all variables and
#'   multiplied by 100.
#' 
#' @param x a matrix of within-sample allele frequencies, as produced by the
#'   function \code{get_wsaf()}.
#' 
#' @return Invisibly returns a list of class `prcomp` with the following 
#'   components
#'   \itemize{
#'       \item{"sdev"}{ the standard deviations of the principal components
#'       (i.e., the square roots of the eigenvalues of the
#'       covariance/correlation matrix, though the calculation is actually done
#'       with the singular values of the data matrix).} \item{"rotation"}{ the
#'       matrix of variable loadings (i.e., a matrix whose columns contain the
#'       eigenvectors). The function \code{princomp} returns this in the element
#'       \code{loadings}.} \item{"center, scale"}{ the centering and scaling
#'       used.} \item{"x"}{ the value of the rotated data (the centred data
#'       multiplied by the rotation matrix). Hence, \code{cov(x)} is the
#'       diagonal matrix \code{diag(sdev^2)}.}
#'       \item{"var"}{ the variance in the data explained by each component.}
#'       \item{"contribution"}{ the percent contribution of a variable (i.e. a
#'       locus) to the overall variation.}
#'       }
#' @importFrom stats prcomp
#' @export

pca_wsaf <- function(x) {
  
  # check inputs
  assert_matrix(x)
  assert_numeric(x)
  if (any(is.na(x))) {
    stop("input matrix cannot contain missing values. If produced using the function get_wsaf() then ensure that imputation is turned on")
  }
  
  # compute PCA
  pca <- prcomp(x)
  
  # compute variance explained
  pca$var <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2) * 100
  
  # compute percent contribution of each locus
  pca$contribution <- abs(pca$rotation)
  pca$contribution <- sweep(pca$contribution, 2, colSums(pca$contribution), "/") * 100
  
  # return
  return(pca)
}

#------------------------------------------------
#' @title Plot variance explained by PCA components
#'
#' @description Plot the variance explained by each PCA component. The number of
#'   components shown is controlled by \code{num_components}, with up to the
#'   first 10 componenets shown by default. If less than the requested number of
#'   components exist, then all the components will be shown.
#'
#' @param pca output of \code{pca_wsaf()} function.
#' @param num_components maximum components to be shown.
#'
#' @importFrom plotly plot_ly layout
#' @export

plot_pca_variance <- function(pca, num_components = 10) {
  
  # potential components
  nc <- min(length(colnames(pca$x)), num_components)
  
  # x has to be factored in order to preserve PC order
  x <- factor(colnames(pca$x)[seq_len(nc)],
              levels = colnames(pca$x)[seq_len(nc)])
  
  # create plot object
  out <- plot_ly(x = x, y = pca$var[seq_len(nc)],
                 type = "bar", width = 500, height = 500) |>
    plotly::layout(title = "Percentage of total variance explaned by PCA",
                   xaxis = list(title = "Principal Component"),
                   yaxis = list(title = "Percentage of variance explained"))
  
  # render and return invisibly
  print(out)
  invisible(out)
}

#------------------------------------------------
#' @title Plot PCA
#'
#' @description Plots either the first 2 or 3 principal components.
#'
#' @param pca output of \code{pca_wsaf()} function.
#' @param num_components numeric for number of components used. Default = 2.
#' @param col vector by which samples are coloured.
#' @param col_palette vector of colours for each group.
#' @param ggplot boolean for plotting using ggplot. Default = FALSE
#'
#' @importFrom plotly plot_ly
#' @importFrom RColorBrewer brewer.pal
#' @export

plot_pca <- function(pca, num_components = 2, col = NULL, col_palette = NULL,
                     ggplot = FALSE) {
  
  # avoid "no visible binding" notes
  PC1 <- PC2 <- NULL
  
  # check inputs
  assert_custom_class(pca, "prcomp")
  assert_in(num_components, c(2,3))
  if (is.null(col)) {
    col <- rep(1, nrow(pca$x))
  }
  assert_vector(col)
  assert_length(col, nrow(pca$x))
  if (is.null(col_palette)) {
    col_palette <- suppressWarnings(brewer.pal(3, "Set1"))
  }
  
  # check num_components
  nc <- min(ncol(pca$x), num_components)
  if (nc != num_components) {
    message("Using all the components (", nc, ") that are available.")
  }
  
  if (num_components == 2) {
    # scatterplot of first 2 principal components
    # 2D scatter
    if (ggplot) {
      
      plot1 <- ggplot2::ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) + 
        geom_point(color = col_palette[col], size = 3) 
      
    } else {
      
      plot1 <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2,
                       color = col, type = "scatter", colors = col_palette,
                       mode = "markers", marker = list(size = 5))
      
    }
  }
  
  if (num_components == 3) {
    # scatterplot of first 3 principal components
    # 3D scatter
    plot1 <- plot_ly(as.data.frame(pca$x[, 1:3]), x = ~PC1, y = ~PC2, z = ~PC3,
                     color = col, type = "scatter3d", colors = col_palette,
                     mode = "markers", marker = list(size = 3))
  }
  
  # render and return invisibly
  print(plot1)
  invisible(plot1)
}

#------------------------------------------------
#' @title Plot PCA contribution of each variable
#'
#' @description Plot PCA contribution of each variable.
#'
#' @param pca output of \code{pca_wsaf()} function.
#' @param component which component to plot.
#' @param chrom the chromosome corresponding to each contribution value.
#' @param pos the genomic position corresponding each contribution value.
#' @param locus_type defines the colour of each bar.
#' @param y_buffer (percent). A buffer added to the bottom of each y-axis,
#'   making room for other annotations to be added.
#'
#' @export

plot_pca_contribution <- function(pca, component = 1, chrom, pos, locus_type = NULL, y_buffer = 0) {
  
  # avoid "no visible binding" notes
  x <- NULL
  
  # check inputs
  assert_custom_class(pca, "prcomp")
  n <- nrow(pca$contribution)
  assert_single_pos_int(component)
  assert_leq(component, n)
  assert_pos_int(chrom)
  assert_length(chrom, n)
  assert_pos_int(pos)
  assert_length(pos, n)
  if (is.null(locus_type)) {
    locus_type <- rep("", n)
  }
  assert_length(locus_type, n)
  assert_single_pos(y_buffer, zero_allowed = TRUE)
  assert_bounded(y_buffer, right = 100)
  
  # extract contributions
  y <- pca$contribution[,component]
  
  # get y ticks and limits
  y_ticks <- pretty(y)
  y_max <- max(y_ticks)
  y_min <- -y_max*y_buffer/100
  
  # load P.falciparum chromosome lengths
  df_chrom_lengths <- Pf_chrom_lengths()
  
  # make dataframes for drawing clustom gridlines
  df_hlines = df_chrom_lengths[rep(1:14, each = length(y_ticks)),]
  df_hlines$y = rep(y_ticks, times = 14)
  
  df_vlines = df_chrom_lengths[rep(1:14, each = 16),]
  df_vlines$x = rep(1:16*2e5 + 1, times = 14)
  df_vlines <- subset(df_vlines, df_vlines$x < df_vlines$length)
  
  # create plotting dataframe
  df <- data.frame(chrom = chrom,
                   pos = pos,
                   locus_type = locus_type,
                   y = y)
  
  # produce basic plot
  plot1 <- ggplot(df) + facet_wrap(~chrom, ncol = 1)
  plot1 <- plot1 + theme(strip.background = element_blank(),
                         strip.text.x = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank())
  
  # add rectangles and grid lines
  plot1 <- plot1 + geom_rect(aes(xmin = 1, xmax = length, ymin = 0, ymax = y_max), col = grey(0.7), size = 0.2, fill = grey(0.95), data = df_chrom_lengths)
  plot1 <- plot1 + geom_segment(aes(x = 1, y = y, xend = length, yend = y), col = grey(0.8), size = 0.1, data = df_hlines)
  plot1 <- plot1 + geom_segment(aes(x = x, y = 0, xend = x, yend = y_max), col = grey(0.8), size = 0.1, data = df_vlines)
  
  # set y scale
  plot1 <- plot1 + scale_y_continuous(breaks = y_ticks, limits = c(y_min,y_max))
  
  # add data
  plot1 <- plot1 + geom_segment(aes(x = pos, y = 0, xend = pos, yend = y, col = locus_type))
  
  # labels and legends
  ylab_title <- paste0("PC", component, " contributions (%)")
  plot1 <- plot1 + xlab("position") + ylab(ylab_title)
  
  # return
  return(plot1)
}

#------------------------------------------------
#' @title PCoA of genomic distances between samples
#'
#' @description Conduct principal coordinate analysis (PCoA) on a matrix of
#'   genomic distances.
#'  
#' @param x matrix of genomic distances, as produced by the function
#'   \code{get_genomic_distance()}.
#'
#' @importFrom ape pcoa
#' @export

pcoa_genomic_distance <- function(x) {
  
  # check inputs
  assert_matrix(x)
  assert_numeric(x)
  
  # mask the diagonal and reflect
  x[is.na(x)] <- 0
  x <- x + t(x)
  
  # compute PCoA
  ret <- ape::pcoa(x)
  
  # return PCoA
  invisible(ret)
}

#------------------------------------------------
#' @title Plot PCoA
#'
#' @description Plots either the first 2 or 3 vectors of PCoA.
#'
#' @param pcoa object of class "pcoa", as produced by \code{pcoa_wsaf()} function.
#' @param num_components numeric for number of components used. Default = 2.
#' @param col vector by which samples are coloured.
#' @param col_palette vector of colours for each group.
#'
#' @importFrom plotly plot_ly
#' @importFrom RColorBrewer brewer.pal
#' @export

plot_pcoa <- function(pcoa, num_components = 2, col = NULL, col_palette = NULL) {
  
  # check inputs
  assert_custom_class(pcoa, "pcoa")
  assert_in(num_components, c(2,3))
  if (is.null(col)) {
    col <- rep(1, nrow(pcoa$vectors))
  }
  assert_vector(col)
  assert_length(col, nrow(pcoa$vectors))
  if (is.null(col_palette)) {
    col_palette <- suppressWarnings(brewer.pal(3, "Set1"))
  }
  
  df <- data.frame(pcoa$vectors[,1:3])
  names(df) <- c("PC1", "PC2", "PC3")
  if (num_components == 2) {
    # scatterplot of first 2 principal components
    # 2D scatter
    plot1 <- plot_ly(df, x = ~PC1, y = ~PC2,
                     color = col, type = "scatter", colors = col_palette,
                     mode = "markers", marker = list(size = 5))
  }
  
  if (num_components == 3) {
    # scatterplot of first 3 principal components
    # 3D scatter
    plot1 <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
                     color = col, type = "scatter3d", colors = col_palette,
                     mode = "markers", marker = list(size = 3))
  }
  
  # render and return invisibly
  print(plot1)
  invisible(plot1)
}
