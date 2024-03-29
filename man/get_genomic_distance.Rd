% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_distance.R
\name{get_genomic_distance}
\alias{get_genomic_distance}
\title{Get genomic distance between samples}
\usage{
get_genomic_distance(x, cutoff = 0.1, report_progress = TRUE)
}
\arguments{
\item{x}{object of class \code{mipanalyzer_biallelic}.}

\item{cutoff}{when calculating weights, correlations below this value are
ignored (see references).}

\item{report_progress}{if \code{TRUE} then a progress bar is printed to the
console.}
}
\description{
Get genomic distance between samples using a distance metric
  that allows for mixed infections and takes account of linkage (see
  references for details).
}
\references{
MalariaGEN Plasmodium falciparum Community Project. "Genomic
  epidemiology of artemisinin resistant malaria". eLIFE (2016).
}
