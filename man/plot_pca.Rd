% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_clustering.R
\name{plot_pca}
\alias{plot_pca}
\title{Plot PCA}
\usage{
plot_pca(
  pca,
  num_components = 2,
  col = NULL,
  col_palette = NULL,
  ggplot = FALSE
)
}
\arguments{
\item{pca}{output of \code{pca_wsaf()} function.}

\item{num_components}{numeric for number of components used. Default = 2.}

\item{col}{vector by which samples are coloured.}

\item{col_palette}{vector of colours for each group.}

\item{ggplot}{boolean for plotting using ggplot. Default = FALSE}
}
\description{
Plots either the first 2 or 3 principal components.
}
