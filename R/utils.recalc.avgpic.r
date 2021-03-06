#' A utility script to recalculate the OneRatioRef, OneRatioSnp, PICRef,	PICSnp, and	AvgPIC by locus after some 
#' populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef,	PICSnp, and	AvgPIC included,
#' but the allelec composition will change when some individuals are removed from the dataset and so the initial statistics will
#' no longer apply. This script recalculates these statistics and places the recalculated values in the appropriate place 
#' in the genlight object.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- dartR:::utils.recalc.avgpic(testset.gl)

utils.recalc.avgpic <- function(x, v=2) {

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n"); stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting utils.recalc.avgpic: Recalculating OneRatioRef, OneRatioSnp, PICRef, PICSnp and AvgPIC\n")
  }
  
  # Do the deed
     t <- as.matrix(x)
     for (i in 1:nLoc(x)) {
       c0 <- length(t[t[,i]==0 & !is.na(t[,i]),i])
       c1 <- length(t[t[,i]==1 & !is.na(t[,i]),i])
       c2 <- length(t[t[,i]==2 & !is.na(t[,i]),i])
       c <- (c0+c1+c2)
       x@other$loc.metrics$OneRatioRef[i] <- (c0+c1)/c
       x@other$loc.metrics$OneRatioSnp[i] <- (c1+c2)/c
       OneRatioRef <- x@other$loc.metrics$OneRatioRef[i]
       OneRatioSnp <- x@other$loc.metrics$OneRatioSnp[i]
       ZeroRatioRef <- 1 - OneRatioRef
       ZeroRatioSnp <- 1 - OneRatioSnp
       x@other$loc.metrics$PICRef[i] <- 1 - ((OneRatioRef*OneRatioRef) + (ZeroRatioRef*ZeroRatioRef))
       x@other$loc.metrics$PICSnp[i] <- 1 - ((OneRatioSnp*OneRatioSnp) + (ZeroRatioSnp*ZeroRatioSnp))
       x@other$loc.metrics$avgPIC[i] <- (x@other$loc.metrics$PICRef[i] + x@other$loc.metrics$PICSnp[i])/2
     }

     if (v > 0) {
       cat("Completed utils.recalc.avgpic\n\n")
     }
   
   return(x)
}