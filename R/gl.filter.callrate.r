#' Filter loci or specimens in a genlight \{adegenet\} object based on call rate
#'
#' SNP datasets generated by DArT have missing values primarily arising from failure to call a SNP because of a mutation
#' at one or both of the the restriction enzyme recognition sites. This script filters out loci (or specimens) for which the call rate is
#' lower than a specified value. The script will also filter out loci (or specimens) in SilicoDArT (presence/absence) datasets where the call rate
#' is lower than the specified value. In this case, the data are missing owing to low coverage.
#'
#' @param gl -- name of the genlight object containing the SNP data, or the genind object containing the SilocoDArT data [required]
#' @param method -- "loc" to specify that loci are to be filtered, "ind" to specify that specimens are to be filtered [default "loc"]
#' @param threshold -- threshold value below which loci will be removed [default 0.95]
#' @return The reduced genlight or genind object, plus a summary
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result2 <- gl.filter.callrate(testset.gl, method="ind", threshold=0.9)



 gl.filter.callrate <- function(gl, method="loc", threshold=0.95) {
 x <- gl
   
  if(class(x) == "genlight") {
     cat("Reporting for a genlight object\n")
   } else if (class(x) == "genind") {
     cat("Reporting for a genind object\n")
   } else {
     cat("Fatal Error: Specify either a genlight or a genind object\n")
     stop()
  }
  cat("Note: Missing values most commonly arise from restriction site mutation.\n\n")
   
  if( method == "loc" ) {
    # Determine starting number of loci and individuals
    n0 <- nLoc(x)
    cat("Initial no. of loci =", n0, "\n")

    if(class(x)=="genlight") {
    # Remove loci with NA count <= 1-threshold
      x2 <- x[ ,glNA(x,alleleAsUnit=FALSE)<=((1-threshold)*nInd(x))]
      x2@other$loc.metrics <- x@other$loc.metrics[glNA(x,alleleAsUnit=FALSE)<=((1-threshold)*nInd(x)),]
      cat ("  No. of loci deleted =", (n0-nLoc(x2)),"\n")

    } else if (class(x)=="genind") {
      x2 <- x[,(colSums(is.na(tab((x))))/nInd(x))<=(1-threshold)]
      idx <- which((colSums(is.na(tab((x))))/nInd(x))<=(1-threshold))
      x2@other$loc.metrics <- x@other$loc.metrics[c(idx),]
      cat ("  No. of loci deleted =", (n0-nLoc(x2)),"\n")
    } else {
      cat("Fatal Error: genlight or genind objects required for call rate filtering!\n"); stop()
    }

  } else if ( method == "ind" ) {
    # Determine starting number of loci and individuals
      n0 <- nInd(x)
      cat("Initial no. of individuals =", n0, "\n")
      # Calculate the individual call rate
      ind.call.rate <- 1 - rowSums(is.na(as.matrix(x)))/nLoc(x)
      # Check that there are some individuals left
      if (sum(ind.call.rate >= threshold) == 0) stop(paste("Maximum individual call rate =",max(ind.call.rate),". Nominated threshold of",threshold,"too stringent.\n No individuals remain.\n"))
      # Extract those individuals with a call rate greater or equal to the threshold
      x2 <- x[ind.call.rate >= threshold,]
      # for some reason that eludes me, this also (appropriately) filters the latlons and the covariates, but see above for locus filtering
      if( class(x) == "genlight") {
        cat ("Filtering a genlight object\n  no. of individuals deleted =", (n0-nInd(x2)), "\nIndividuals retained =", nInd(x2),"\n")
      }
      if( class(x) == "genind") {
        cat ("Filtering a genind object\n  No. of individuals deleted =", (n0-nInd(x2)), "\n  Individuals retained =", nInd(x2),"\n")
      }
      # Report individuals that are excluded on call rate
      if (any(ind.call.rate <= threshold)) {
        x3 <- x[ind.call.rate <= threshold,]
        if (length(x3) > 0) {
          cat("List of individuals deleted because of low call rate\n",indNames(x3),"\n")
          cat("   from populations\n",as.character(pop(x3)),"\n")
        }
      }  
  }  else {
      cat("Fatal Error: the method parameter must be specified as either loc or ind\n")
      stop()
  }

   # REPORT A SUMMARY
   cat("Summary of filtered dataset\n")
   cat(paste("  Call Rate >",threshold,"\n"))
   cat(paste("  No. of loci:",nLoc(x2),"\n"))
   cat(paste("  No. of individuals:", nInd(x2),"\n"))
   cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))

    return(x2)
}

 