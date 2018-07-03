#' Identify loci that are sex linked in specimens in a genlight \{adegenet\} object
#'
#' Alleles unique to the Y or W chromosome and monomorphic on the X chromosomes will appear in the SNP dataset as 
#' genotypes that are heterozygotic in all individuals of the heterogametic sex and homozygous in all individuals 
#' of the homogametic sex.
#' 
#' This script will identify loci with alleles that behave in this way, as putative sex specific SNP markers.
#' 
#' Sex of the individuals for which sex is known with certainty is to be held in the variable x@other$ind.metrics$sex, 
#' as M for male, F for female, NA otherwise. The script abbreviates the entries here to the first character. So coding of "Female" and "Male" works as well. Character are also converted to upper cases.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param t.het -- tolerance, that is t.het = 0.05 means that 5% of the homogametic sex can be heterozygous and still 
#' be regarded as consistent with a sex specific marker [default 0]
#' @param t.hom -- tolerance, that is t.hom = 0.05 means that 5% of the heterogametic sex can be homozygous and still
#' be regarded as consistent with a sex specific marker [default 0]
#' @param verbose -- verbosity: logical indicating whether outputs should be printed to the console (default: FALSE)
#' @param na.rm -- logical: should NAs in sex assignments be ignored?
#' @return The list of sex specific loci
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- gl.sexlinkage(testset.gl)

gl.sexlinkage <- function(x, t.het = 0, t.hom = 0, verbose = FALSE, na.rm = FALSE) {

  # remove NAs from data set
  if (na.rm & any(is.na(x@other$ind.metrics$sex))) {
    x <- x[!is.na(x@other$ind.metrics$sex)]
  }
  
  if (verbose) {
    cat("Starting gl.sexlinkage: Identifying sex linked loci\n")
  }
  
  if(class(x) != "genlight") {
    cat("Fatal Error: genlight object required for gl.sexlinkage!\n"); stop("Execution terminated\n")
  }

  # Extract the sex variable from whereever it may be -- might need alteration    
  sex <- toupper(substr(x@other$ind.metrics$sex, 1, 1))

  # Extract the data for the females
  matf <- as.matrix(x[sex == "F"])
  # For each individual
    f <- array(data = NA, dim = c(ncol(matf), 3))
    for (i in seq_len(ncol(matf))) {
      for (j in 1:3) {
        f[i, j] <- length(which(matf[, i] == (j - 1)))
      }  
    }
    dff <- data.frame(f)
    row.names(dff) <- locNames(x)
    colnames(dff) <- c("F0", "F1", "F2")

  # Extract the data for the males
  matm <- as.matrix(x[sex == "M"])
  
  # For each individual
  m <- array(data = NA, dim = c(ncol(matm), 3))
  for (i in seq_len(ncol(matm))) {
    for (j in 1:3) {
      m[i, j] <- length(which(matm[, i] == (j - 1)))
    }  
  }
  dfm <- data.frame(m)
  row.names(dfm) <- locNames(x)
  colnames(dfm) <- c("M0", "M1", "M2")

  # Combine the two files
  Trimmed_Sequence <- x@other$loc.metrics$TrimmedSequence
  df <- cbind(dff, dfm, Trimmed_Sequence)
  a <- strsplit(row.names(df), split = "-")
  a <- do.call(rbind, a)
  a <- as.numeric(a[, 2])
  
  df$Trimmed_Sequence <- as.character(df$Trimmed_Sequence)
  b <- substr(df$Trimmed_Sequence, 1, a)
  c <- substr(df$Trimmed_Sequence, a + 1, a + 1)
  c <- tolower(c)
  d <- substr(df$Trimmed_Sequence, a + 2, nchar(df$Trimmed_Sequence))
  
  df$Trimmed_Sequence <- paste0(b, c, d)
  df$AvgCountRef <- x@other$loc.metrics$AvgCountRef
  df$AvgCountSnp <- x@other$loc.metrics$AvgCountSnp
  
  # Check for hets in all males, homs in all females (XY); ditto for ZW
  sumf <- df$F0 + df$F1 + df$F2
  summ <- df$M0 + df$M1 + df$M2
  df_zero_sum <- df[(sumf == 0) | (summ == 0), ]
  df <- df[(sumf > 0) & (summ > 0), ]

  # Locus present on both Z and W that differentiates males and females
  #   (only hets in females, only homs in males)
  zw <- df[(df$F1 / (sumf[(sumf > 0) & (summ > 0)])) >= (1 - t.hom) &
             (df$M1 / (summ[(sumf > 0) & (summ > 0)])) <= (0 + t.het), ]
  
  # Locus present on both X and Y that differentiates males and females
  #    (only homs in females, only hets in males)
  xy <- df[(df$F1 / (sumf[(sumf > 0) & (summ > 0)])) <= (0 + t.het) &
             (df$M1 / (summ[(sumf > 0) & (summ > 0)])) >= (1 - t.hom), ]

  # Z-linked/only on Z (only homs in females, hets in males)
  zl <- df[(df$F1 / sumf[(sumf > 0) & (summ > 0)]) <= (0 + t.het) &
             (df$M1 / (summ[(sumf > 0) & (summ > 0)])) >= (1 - t.hom), ]
  
  # W-linked/only on W (only homs in females, absent in males)
  wl <- df_zero_sum[(df_zero_sum$F1 / sumf[(sumf == 0) | (summ == 0)]) <= (0 + t.het) &
                      summ[(sumf == 0) | (summ == 0)] == 0, ]
  
  # X-linked/only on X (hets in females, only homs in males)
  xl <- df[(df$M1 / summ[(sumf > 0) & (summ > 0)]) <= (0 + t.het) &
             (df$F1 / sumf[(sumf > 0) & (summ > 0)]) >= (1 - t.hom), ]
  
  # Y-linked/only on Y (absent in females, only homs in males)
  yl <- df_zero_sum[(df_zero_sum$M1 / summ[(sumf == 0) | (summ == 0)]) <= (0 + t.het) &
                      sumf[(sumf == 0) | (summ == 0)] == 0, ]
  
  # print only if verbose
  if (verbose) {
    
    # Locus present on both Z and W that differentiates males and females
    #   (only hets in females, only homs in males)
    if (nrow(zw) == 0) {
      cat("No loci present on Z and W that differentiate males and females (ZZ/ZW)\n")
    } else {
      cat("\nFound loci on Z and W that differentiate males and females (ZZ/ZW)\n")
      cat(paste("  Threshold proportion for homozygotes in the heterozygotic sex (ZW)", t.hom, ";\n")) 
      cat(paste("  for heterozygotes in the homozygotic sex (ZZ)", t.het, "\n"))
      cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
      print(zw)
      cat("Note: SNP location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
      cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
    }
    
    # Locus present on both X and Y that differentiates males and females
    #    (only homs in females, only hets in males)
    if (nrow(xy) == 0) {
      cat("No loci present on X and Y that differentiate males and females (XX/XY)\n")
    } else {
      cat("\nFound loci on X and Y that differentiate males and females (XX/XY)\n")
      cat(paste("  Threshold proportion for homozygotes in the heterozygotic sex (XY)",t.hom,";\n")) 
      cat(paste("  for heterozygotes in the homozygotic sex (XX)",t.het,"\n"))
      cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
      print(xy)
      cat("Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
      cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
    }
    
    # Z-linked/only on Z (only homs in females, hets in males)
    if (nrow(zl) == 0) {
      cat("No loci present that are Z-linked or found only on Z")
    } else {
      cat("\nFound loci that are Z-linked or found only on Z\n")
      cat(paste("  Threshold proportion for heterozygotes in the homozygotic sex (W)", t.het, ";\n")) 
      cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
      print(zl)
      cat("Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
      cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
    }
    
    # W-linked/only on W (only homs in females, absent in males)
    if (nrow(wl) == 0) {
      cat("No loci present that are W-linked or found only on W")
    } else {
      cat("\nFound loci that are W-linked or found only on W\n")
      cat(paste("  Threshold proportion for heterozygotes in the homozygotic sex (W)", t.het, ";\n")) 
      cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
      print(wl)
      cat("Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
      cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
    }
    
    # X-linked/only on X (hets in females, only homs in males)
    if (nrow(xl) == 0) {
      cat("No loci present that are X-linked or found only on X")
    } else {
      cat("\nFound loci that are X-linked or found only on X\n")
      cat(paste("  Threshold proportion for heterozygotes in the homozygotic sex (Y)", t.het, ";\n")) 
      cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
      print(xl)
      cat("Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
      cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
    }
    
    # Y-linked/only on Y (absent in females, only homs in males)
    if (nrow(yl) == 0) {
      cat("No loci present that are Y-linked or found only on Y")
    } else {
      cat("\nFound loci that are Y-linked or found only on Y\n")
      cat(paste("  Threshold proportion for heterozygotes in the homozygotic sex (Y)", t.het, ";\n")) 
      cat("  0 = homozygous reference; 1 = heterozygous; 2 = homozygous alternate\n")
      print(yl)
      cat("Note: Snp location in Trimmed Sequence indexed from 0 not 1, SNP position in lower case\n")
      cat("Note: The most reliable putative markers will have AvgCount for Ref or Snp 10 or more, one ca half the other\n")
    }
    
  }
  
  out <- list(zw, xy, zl, wl, xl, yl)
  
  class(out) <- 'sexlinkage'
  
  out
  
}

print.sexlinkage <- function(x, ...) {
  
  UseMethod(print, x)
  
}
 
print.sexlinkage.default <- function(x) {
  
  NULL
  
}