#' @rdname sexlinkage
#' 
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
#' @param x
#' @param object
#' 
#' @return The list of sex specific loci
#' 
#' @export
#' 
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' 
#' @examples
#' result <- gl.sexlinkage(testset.gl)

gl.sexlinkage <- function(x,
                          t.het = 0, t.hom = 0,
                          verbose = FALSE,
                          na.rm = FALSE) {

  # remove NAs from data set
  if (na.rm & any(is.na(x@other$ind.metrics$sex))) {
    x <- x[!is.na(x@other$ind.metrics$sex)]
  }
  
  if (verbose) {
    cat("Starting gl.sexlinkage: Identifying sex linked loci\n")
  }
  
  if(class(x) != "genlight") {
    stop("Fatal Error: genlight object required for gl.sexlinkage\n",
         call. = FALSE)
  }

  # Extract the sex variable from whereever it may be -- might need alteration    
  sex <- toupper(substr(x@other$ind.metrics$sex, 1, 1))

  # Extract the data for the females
  matf <- as.matrix(x[sex == "F"])
  
  # For each individual
  f <- matrix(NA, nrow = ncol(matf), ncol = 3)
  for (i in 1:3) {
    f[, i] <- apply(matf, 2, function(x) sum(x == (i - 1), na.rm = TRUE))
  }
  dff <- data.frame(f)
  row.names(dff) <- locNames(x)
  colnames(dff) <- c("F0", "F1", "F2")
  
  # Extract the data for the males
  matm <- as.matrix(x[sex == "M"])
  
  # For each individual
  m <- matrix(NA, nrow = ncol(matm), ncol = 3)
  for (i in 1:3) {
    m[, i] <- apply(matm, 2, function(x) sum(x == (i - 1), na.rm = TRUE))
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

  # W-linked/only on W (only homs in females, absent in males)
  wl <- df_zero_sum[(df_zero_sum$F1 / sumf[(sumf == 0) | (summ == 0)]) <= (0 + t.het) &
                      summ[(sumf == 0) | (summ == 0)] == 0, ]
  
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
  
  out <- list(fem_het_male_hom = zw,
              fem_hom_male_het = xy,
              fem_hom_male_absent = wl,
              fem_absent_male_hom = yl)

  out <- as.sexlinkage(out)
  
  out
  
}


#' @rdname sexlinkage
#'
#' @export
#' 
#' @examples
#'
#' # check if an object is a sexlinkage object
#'   
#' \dontrun{
#' is.sexlinkage(object)
#' }

is.sexlinkage <- function (object) {
  inherits(object, 'sexlinkage')
}

#' @rdname sexlinkage
#'
#' @export
#'
#' @examples
#' 
#' # Print information about a 'sexlinkage' object
#'
#' \dontrun{
#' print(object)
#' }

print.sexlinkage <- function (x, ..., system = c('xy', 'zw')) {
  
  # work out if any loci are sex-linked
  num_loci <- c(nrow(x$fem_het_male_hom),
                nrow(x$fem_hom_male_het),
                nrow(x$fem_hom_male_absent),
                nrow(x$fem_absent_male_hom))
  
  if (any(num_loci) > 0) {
    if (is.null(system) | !(system %in% c('xy', 'zw'))) {
      if (nrow(x$fem_het_male_hom)) {    
        cat(paste0('Identified ', nrow(x$fem_het_male_hom),
                   ' loci with evidence of heterogametic something',
                   '\n'))
      }
      if (nrow(x$fem_hom_male_het)) {    
        cat(paste0('Identified ', nrow(x$fem_hom_male_het),
                   ' loci with evidence of male hemizygosity',
                   '\n'))
      }
      if (nrow(x$fem_hom_male_absent)) {    
        cat(paste0('Identified ', nrow(x$fem_hom_male_absent),
                   ' loci with evidence of male hemizygosity',
                   '\n'))
      }
      if (nrow(x$fem_absent_male_hom)) {    
        cat(paste0('Identified ', nrow(x$fem_absent_male_hom),
                   ' loci with evidence of male hemizygosity',
                   '\n'))
      }
    }
    if (system %in% c('xy', 'zw')) {
      if (system == 'xy') { 
        if (nrow(x$fem_het_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_het_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
        if (nrow(x$fem_hom_male_het)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_het),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
        if (nrow(x$fem_hom_male_absent)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_absent),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
        if (nrow(x$fem_absent_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_absent_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
      }
      if (system == 'zw') {
        if (nrow(x$fem_het_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_het_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
        if (nrow(x$fem_hom_male_het)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_het),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
        if (nrow(x$fem_hom_male_absent)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_absent),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
        if (nrow(x$fem_absent_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_absent_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
        }
      }
    }
  } else {
    cat(paste0('No sex-linked loci were identified\n'))
  } 
}

#' @rdname sexlinkage
#'
#' @export
#'
#' @examples
#' 
#' # Summarise a 'sexlinkage' object
#'
#' \dontrun{
#' summary(object)
#' }

summary.sexlinkage <- function (object, ..., system = c('xy', 'zw')) {
  
  # work out if any loci are sex-linked
  num_loci <- c(nrow(x$fem_het_male_hom),
                nrow(x$fem_hom_male_het),
                nrow(x$fem_hom_male_absent),
                nrow(x$fem_absent_male_hom))
  
  if (any(num_loci) > 0) {
    if (is.null(system) | !(system %in% c('xy', 'zw'))) {
      if (nrow(x$fem_het_male_hom)) {    
        cat(paste0('Identified ', nrow(x$fem_het_male_hom),
                   ' loci with evidence of heterogametic something',
                   '\n'))
        print(x$fem_het_male_hom)
      }
      if (nrow(x$fem_hom_male_het)) {    
        cat(paste0('Identified ', nrow(x$fem_hom_male_het),
                   ' loci with evidence of male hemizygosity',
                   '\n'))
        print(x$fem_hom_male_het)
      }
      if (nrow(x$fem_hom_male_absent)) {    
        cat(paste0('Identified ', nrow(x$fem_hom_male_absent),
                   ' loci with evidence of male hemizygosity',
                   '\n'))
        print(x$fem_hom_male_absent)
      }
      if (nrow(x$fem_absent_male_hom)) {    
        cat(paste0('Identified ', nrow(x$fem_absent_male_hom),
                   ' loci with evidence of male hemizygosity',
                   '\n'))
        print(x$fem_absent_male_hom)
      }
    }
    if (system %in% c('xy', 'zw')) {
      if (system == 'xy') { 
        if (nrow(x$fem_het_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_het_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_het_male_hom)
        }
        if (nrow(x$fem_hom_male_het)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_het),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_hom_male_het)
        }
        if (nrow(x$fem_hom_male_absent)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_absent),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_hom_male_absent)
        }
        if (nrow(x$fem_absent_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_absent_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_absent_male_hom)
        }
      }
      if (system == 'zw') {
        if (nrow(x$fem_het_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_het_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_het_male_hom)
        }
        if (nrow(x$fem_hom_male_het)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_het),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_hom_male_het)
        }
        if (nrow(x$fem_hom_male_absent)) {    
          cat(paste0('Identified ', nrow(x$fem_hom_male_absent),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_hom_male_absent)
        }
        if (nrow(x$fem_absent_male_hom)) {    
          cat(paste0('Identified ', nrow(x$fem_absent_male_hom),
                     ' loci with evidence of male hemizygosity',
                     '\n'))
          print(x$fem_absent_male_hom)
        }
      }
    }
  } else {
    cat(paste0('No sex-linked loci were identified\n'))
  } 
  
}

# internal function: create sexlinkage object
as.sexlinkage <- function (object) {
  as_class(object, name = 'sexlinkage', type = 'list')
}

# set an object class
as_class <- function (object, name, type = c("function", "list")) {
  
  type <- match.arg(type)
  stopifnot(inherits(object, type))
  class(object) <- c(name, class(object))
  
  object
  
}
