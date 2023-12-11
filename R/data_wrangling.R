#------------------------------------------------
#' @title Load system file
#'
#' @description Load a file from within the inst/extdata folder of the
#'   MIPanalyzer package. File extension must be one of .csv, .txt, or .rds.
#'
#' @param name the name of a file within the inst/extdata folder.
#'
#' @importFrom data.table fread
#' @export

mipanalyzer_file <- function(name) {
  
  # check that valid file extension
  ext <- strsplit(name, "\\.")[[1]]
  ext <- ext[length(ext)]
  assert_in(ext, c("txt", "csv", "rds"), message = "file extension not valid")
  
  # get full file path
  name_full <- system.file("extdata/", name, package='MIPanalyzer', mustWork = TRUE)
  
  # read in file
  if (ext == "rds") {
    ret <- readRDS(name_full)
  } else {
    ret <- data.table::fread(name_full, data.table = FALSE)
  }
  
  return(ret)
}

#------------------------------------------------
#' @title Convert vcf to biallelic mipanalyzer data class
#'
#' @description Convert vcf to biallelic mipanalyzer data class.
#'
#' @param file path to vcf file.
#' @param vcfR object of class \code{vcfR}.
#' @param verbose if reading from file, whether to read in verbose manner.
#'
#' @export

vcf_to_mipanalyzer_biallelic <- function(file = NULL, vcfR = NULL, verbose = TRUE) {
  
  # avoid "no visible binding" notes
  CHROM <- POS <- QUAL <- n <- locicount <- NULL
  
  # check inputs
  if (!xor(!is.null(file), !is.null(vcfR))) {
    stop("Must specify one input: either a raw vcf file path or a vcfR object")
  }
  
  # get vcf object
  if (!is.null(vcfR)){
    assert_custom_class(vcfR, "vcfR")
    vcf <- vcfR
  } else {
    assert_file_exists(file)
    vcf <- vcfR::read.vcfR(file = file, verbose = verbose)
  }
  
  message("Processing")
  
  # check that vcf is normalised
  vcf_unnormalised <- vcf@fix |>
    tibble::as_tibble() |>
    dplyr::group_by(CHROM, POS) |>
    dplyr::summarise(locicount = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::pull(locicount) |>
    max() != 1
  if (vcf_unnormalised) {
    stop("This is not a normalized vcf. Consider running bcftools norm, and/or review how the vcf was created")
  }
  
  # check that biallelic at all loci
  if (!all(vcfR::is.biallelic(vcf))) {
    stop("All loci must be biallelic")
  }
  
  # extract coverage and counts matrices
  coverage <- t(vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE))
  counts_raw <- t(vcfR::extract.gt(vcf, element = "AD"))
  counts <- vcfR::masplit(counts_raw, record = 1, sort = FALSE, decreasing = FALSE)
  
  # check that all missing fields correspond between coverage and counts
  if (!identical(is.na(coverage), is.na(counts))) {
    stop("Cannot have fields that are NA in coverage but not in counts, or vice versa")
  }
  
  # extract sample info
  samples <- data.frame(SAMPLE_ID = colnames(vcf@gt)[-1], stringsAsFactors = FALSE)
  
  # extract loci and specify some columns classes
  loci <- vcf@fix |>
    as.data.frame() |>
    dplyr::mutate(POS = as.numeric(as.character(POS)),
                  QUAL = as.numeric(as.character(QUAL)))
  
  # keep vcf meta for downstream processes
  meta <- vcf@meta
  
  # initialise filter history
  filter_history <- data.frame(description = "raw data",
                               samples = nrow(coverage),
                               loci = ncol(coverage),
                               n_missing = sum(is.na(coverage)),
                               prop_missing = mean(is.na(coverage)),
                               function_call = NA,
                               stringsAsFactors = FALSE)
  
  # create return list
  ret <- list(coverage = coverage,
              counts = counts,
              samples = samples,
              loci = loci,
              filter_history = filter_history,
              vcfmeta = meta)
  
  message("Done")
  
  # return in mipanalyzer_biallelic class
  class(ret) <- "mipanalyzer_biallelic"
  return(ret)
}

#------------------------------------------------
#' @title Convert vcf to multiallelic mipanalyzer data class
#'
#' @description Convert vcf to multiallelic mipanalyzer data class.
#'
#' @param file path to vcf file.
#' @param vcfR object of class \code{vcfR}.
#' @param verbose if reading from file, whether to read in verbose manner.
#' @export

vcf_to_mipanalyzer_multiallelic <- function(file = NULL, vcfR = NULL, verbose = TRUE) {
  
  # avoid "no visible binding" notes
  CHROM <- POS <- QUAL <- n <- locicount <- NULL
  
  # check inputs
  if (!xor(!is.null(file), !is.null(vcfR))) {
    stop("Must specify one input: either a raw vcf file path or a vcfR object")
  }
  
  # get vcf object
  if (!is.null(vcfR)){
    assert_custom_class(vcfR, "vcfR")
    vcf <- vcfR
  } else {
    assert_file_exists(file)
    vcf <- vcfR::read.vcfR(file = file, verbose = verbose)
  }
  
  message("Processing")
  
  # check that vcf is normalised
  vcf_unnormalised <- vcf@fix |>
    tibble::as_tibble() |>
    dplyr::group_by(CHROM, POS) |>
    dplyr::summarise(locicount = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::pull(locicount) |>
    max() != 1
  if (vcf_unnormalised) {
    stop("This is not a normalized vcf. Consider running bcftools norm, and/or review how the vcf was created")
  }
  
  # extract counts into array and compute coverage as sum of counts
  counts_raw <- t(vcfR::extract.gt(vcf, element = "AD"))
  counts <- array(NA, dim = c(4, nrow(counts_raw), ncol(counts_raw)))
  for (i in 1:4) {
    counts[i,,] <- vcfR::masplit(counts_raw, record = i, sort = FALSE)
  }
  coverage <- colSums(counts, na.rm = TRUE)
  coverage[coverage == 0] <- NA
  
  # extract sample info
  samples <- data.frame(SAMPLE_ID = colnames(vcf@gt)[-1], stringsAsFactors = FALSE)
  
  # extract loci and specify some columns classes
  loci <- vcf@fix |>
    as.data.frame() |>
    dplyr::mutate(POS = as.numeric(as.character(POS)),
                  QUAL = as.numeric(as.character(QUAL)))
  
  # keep vcf meta for downstream processes
  meta <- vcf@meta
  
  # initialise filter history
  filter_history <- data.frame(description = "raw data",
                               samples = nrow(coverage),
                               loci = ncol(coverage),
                               n_missing = sum(is.na(coverage)),
                               prop_missing = mean(is.na(coverage)),
                               function_call = NA,
                               stringsAsFactors = FALSE)
  
  # create return list
  ret <- list(coverage = coverage,
              counts = counts,
              samples = samples,
              loci = loci,
              filter_history = filter_history,
              vcfmeta = meta)
  
  message("Done")
  
  # return in mipanalyzer_biallelic class
  class(ret) <- "mipanalyzer_multiallelic"
  return(ret)
}

#------------------------------------------------
#' @title Converts a MIPanalyzer biallelic object to vcfR
#'
#' @description Converts an object of class \code{mipanalyzer_biallelic} to
#'   \code{vcfR} format.
#'
#' @param input an object of class \code{mipanalyzer_biallelic}.
#' @param cutoff the within-sample non-referent allele frequency cutoff to
#'   transform your biallelic site to a genotype matrix.
#'
#' @export

mipanalyzer_biallelic_to_vcfR <- function(input = NULL, cutoff = 0.1) {
  
  # avoid "no visible binding" notes
  CHROM <- POS <- ID <- REF <- ALT <- QUAL <- FILTER <- INFO <- NULL
  
  if(!inherits(input, c("mipanalyzer_biallelic"))) {
    stop("This function only works on objects of class mipanalyzer_biallelic")
  }
  
  # setup for gt
  wsraf <- input$counts / input$coverage
  GT <- ifelse(wsraf > 1 - cutoff, "0/0",
               ifelse(wsraf < cutoff, "1/1",
                      ifelse(!is.na(wsraf), "0/1", NA)))
  ADref <- input$counts
  ADalt <- input$coverage - input$counts
  DP <- input$coverage
  
  gt <- t(ifelse(is.na(GT), NA, sprintf("%s:%s,%s:%s", GT, ADref, ADalt, DP)))
  
  # append format column and sample names
  gt <- cbind(FORMAT = "GT:AD:DP", gt)
  colnames(gt)[2:ncol(gt)] <- input$samples$SAMPLE_ID
  
  # getFix
  fix <- input$loci |>
    dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) |>
    as.matrix()
  fix[,"POS"] <- stringr::str_replace_all(fix[,"POS"], "\\s", "")  
  
  # get meta
  meta <- append(input$vcfmeta, "##MIPanalyzer=This vcf was filtered and modified by the MIPanalyzer R package")
  
  # write out new vcfRobj
  ret <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(ret)
}
