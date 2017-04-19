#' Test to see if two populations are fixed at a given locus
#'
#' This script compares two percent allele frequencies
#' and reports TRUE if they represent a fixed difference, FALSE otherwise.
#'
#' A fixed difference at a locus occurs when two populations share no alleles, noting that SNPs are biallelic (ploidy=2).
#' Tollerance in the definition of a fixed difference is provided by the t parameter. For example, t=0.05 means
#' that SNP allele frequencies of 95,5 and 5,95 percent will be reported as fixed (TRUE).
#'
#' @param s1 -- percentage SNP allele frequency for the first population [required]
#' @param s2 -- percentage SNP allele frequency for the second population [required]
#' @param t -- threshold value for tollerance in when a difference is regarded as fixed [default 0]
#' @export
#' @return TRUE (fixed difference) or FALSE (alleles shared) or NA (one or both s1 or s2 missing)
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' 
#' is.fixed(0,100)
#' is.fixed(100,0)
#' is.fixed(80,0)
#' is.fixed(100,NA)
#' is.fixed(0,NA)
#' is.fixed(NA,0)
#' is.fixed(NaN,100)
#' is.fixed(NaN,0)
#' is.fixed(100,NaN)
#' is.fixed(0,NaN)
#' is.fixed(NaN,NaN)
#' is.fixed(NA,NA)
#' is.fixed(s1=100, s2=0, t=0)
#' is.fixed(96, 4, t=0.05)
#' 
#' @seealso \code{\link{gl.fixed.diff}}

is.fixed<-function(s1, s2, t=0){

  if (is.na(s1) | is.na(s2)) {
    result<-NA
  } else {
    result<-0
    if ((s1<=t*100) & (s2>=(100-t*100))) {result<-1}
    if ((s2<=t*100) & (s1>=(100-t*100))) {result<-1}
  }
  return(result)
}


