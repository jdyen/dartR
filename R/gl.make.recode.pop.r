#' Create a proforma recode_ind_table file for reassigning individual (=specimen) names
#'
#' Renaming individuals may be required when there have been errors in labelling arising
#' in the process from sample to DArT files. There may be occasions where renaming
#' individuals is required for preparation of figures. Caution needs to be exercised
#' because of the potential for breaking the "chain of evidence" between the samples themselves
#' and the analyses. REcoding individuals can be done with a recode table (csv).
#'
#' This script facilitates 
#' the construction of a recode table by producing a proforma file with
#' current individual (=specimen) names in two identical columns. Edit the second
#' column to reassign individual names. Use keyword Delete to delete an individual.
#' 
#' Apply the recoding using gl.recode.ind(). Deleting individuals
#' can potentially generate monomorphic loci or loci with all
#' values missing. Clean this up with gl.filter.monomorphic().
#'
#' @param x -- name of the genlight object containing the SNP data, or the genind object containing the SilocoDArT data [required]
#' @param outfile -- name of the new proforma file [default recode_ind_table.csv]
#' @return A vector containing the new individual names
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' result <- gl.make.recode.ind(gl, outfile="Emmac_recode_ind.csv")
#' }

 gl.make.recode.pop <- function(x, outfile="recode_ind_table.csv") {

 mat <- cbind(levels(pop(x)),levels(pop(x)))
 write.table(mat, file=outfile, sep=",", row.names=FALSE, col.names=FALSE)
 
 cat(paste("Proforma recode table written to:",outfile,"\n"))
 
 return(levels(pop(x)))
 
 }