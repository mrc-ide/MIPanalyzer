# Figured you would want to tweak this slightly before putting them into main (or a different "out" file)?
# Also open to discussion on the assign GT. That is just someting quick I use for the RMCL .... but probably should be more robust...




#' @title Assign the GT calls from the REFERENT WSAF
#' @param wsraf A within-sample non-referent allele frequncy from biallelic SNPs as a \code{matrix}.  
#' @param cutoff The allele-frequency cutoff to determine genotype calls
#' @description This function takes in the within-sample non-referent allele frequency matrix and converts it to a (diploid) genotype matrix. The genotype matrix is 
#' character matrix with "0/0", "0/1", "1/1" representing the homozygote referent, heterozygote, and homozygote alternative calls. The genotype call is 
#' determined by the allele frequency and the \code{cutoff} set by the user. Specifically, an allele frequency less than the cutoff or greater than \code{1-cutoff} correspond to a
#' homozygote alternative and homozygote referent call, respectively. All other allele-frequencies will be converted to the heterozygote call \code{"0/1"}. 
#' 
#' @return A GT matrix as a \code{matrix}.
#' not exported

assignGTfrombiWSRAF <- function(wsraf, cutoff = 0.1){
  
  # extract GT matrix
  GT <- matrix(NA, dim(wsraf)[1], dim(wsraf)[2])
  
  # round accordin to cutoff
  GT <- ifelse(wsraf > 1-cutoff, "0/0",
               ifelse(wsraf < 0+cutoff, "1/1",
                      ifelse(!is.na(wsraf), "0/1", NA)))
  
  return(GT)
}


#' @title Convert a MIPanalyzer Biallelic Object to a vcfR object
#' @param input is of class \code{mipanalyzer_biallelic}
#' @param cutoff is the within-sample non-referent allele frequency cutoff to transform your biallelic site to a genotype matrix.
#'
#' @export


MIPanalyzerbi2vcfR <- function(input = NULL, cutoff = 0.1){
  
  if(!inherits(input, c("mipanalyzer_biallelic"))){
    stop("This function only works on objects of class mipanalyzer_biallelic or mipanalyzer_multiallelic, not class ", class(mipobj))
  }
  
  # setup for gt
  wsraf <- input$counts/input$coverage
  GT <- assignGTfrombiWSRAF(wsraf, cutoff = cutoff)
  ADref <- input$counts
  ADalt <- input$coverage - input$counts
  DP <- input$coverage
  
  gt <- t( ifelse(is.na(GT), ".", paste0(GT, ":", ADref, ",", ADalt, ":", DP)) ) # loci as rows, smpls as columns
  
  # append format column and sample names
  gt <- cbind(FORMAT = "GT:AD:DP", gt)
  colnames(gt)[2:ncol(gt)] <- input$samples$ID
  
  # getFix
  fix <- as.matrix(input$loci[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]) # must be in this order and only these
  fix[,2] <- stringr::str_replace_all(fix[,2], "\\s", "")  
  
  # get meta
  meta <- append(input$vcfmeta, "##MIPanalyzer=This vcf was filtered and modified by the MIPanalyzer R package")
  
  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)
  
}

