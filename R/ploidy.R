#' Estimate baseline ploidy from allele frequencies
#'
#' Estimates population-level ploidy from pseudo-bulk allele-frequency
#' distributions. Triploid populations are identified from enrichment
#' around AF = 1/3 and 2/3. Tetraploid populations require evidence
#' around AF = 1/4 and 3/4.
#'
#' @param TapestriExperiment TapestriExperiment object.
#' @param sample.feature colData column defining cell populations.
#' @param min.cells Minimum cells with a measurement for a variant.
#' @param min.variants Minimum informative variants required.
#' @param tolerance Distance from theoretical allele fraction.
#' @param known.ploidy Optional named vector of known ploidies.
#'
#' @return data.frame with ploidy inference metrics.
#' @export

estimatePloidyFromAF <- function(
    TapestriExperiment,
    sample.feature = "cluster",
    min.cells = 20,
    min.variants = 20,
    tolerance = 0.06,
    known.ploidy = NULL
) {
  
  if (!"alleleFrequency" %in%
      SingleCellExperiment::altExpNames(TapestriExperiment)) {
    
    cli::cli_abort(
      "alleleFrequency altExp not found."
    )
  }
  
  af.exp <- SingleCellExperiment::altExp(
    TapestriExperiment,
    "alleleFrequency"
  )
  
  af <- SummarizedExperiment::assay(
    af.exp,
    "alleleFrequency"
  )
  
  # Mission Bio may store AF as 0-100
  if (max(af, na.rm = TRUE) > 1.5) {
    af <- af / 100
  }
  
  groups <- as.character(
    SummarizedExperiment::colData(
      TapestriExperiment
    )[[sample.feature]]
  )
  
  names(groups) <- colnames(TapestriExperiment)
  
  groups <- groups[colnames(af)]
  
  group.names <- unique(groups)
  
  results <- lapply(
    group.names,
    function(current.group) {
      
      cells <- which(
        groups == current.group
      )
      
      x <- af[
        ,
        cells,
        drop = FALSE
      ]
      
      n.obs <- rowSums(
        !is.na(x)
      )
      
      pseudo.bulk.af <-
        matrixStats::rowMedians(
          x,
          na.rm = TRUE
        )
      
      keep <- (
        n.obs >= min.cells &
          pseudo.bulk.af >= 0.10 &
          pseudo.bulk.af <= 0.90
      )
      
      vaf <- pseudo.bulk.af[keep]
      
      n.variants <- length(vaf)
      
      if (n.variants < min.variants) {
        
        return(
          data.frame(
            population = current.group,
            n.variants = n.variants,
            support.2n = NA,
            support.3n = NA,
            support.4n.outer = NA,
            estimated.ploidy = NA,
            call = "low_information",
            stringsAsFactors = FALSE
          )
        )
      }
      
      
      # ======================================================
      # DISTANCE FROM THEORETICAL VAF STATES
      # ======================================================
      
      d2 <- abs(
        vaf - 1 / 2
      )
      
      d3 <- pmin(
        abs(vaf - 1 / 3),
        abs(vaf - 2 / 3)
      )
      
      # Evidence specific for tetraploidy.
      #
      # 0.5 cannot distinguish CN2 AB from CN4 AABB.
      # Therefore only 0.25 / 0.75 are considered
      # tetraploid-specific evidence.
      d4.outer <- pmin(
        abs(vaf - 1 / 4),
        abs(vaf - 3 / 4)
      )
      
      
      support.2n <- mean(
        d2 <= tolerance
      )
      
      support.3n <- mean(
        d3 <= tolerance
      )
      
      support.4n.outer <- mean(
        d4.outer <= tolerance
      )
      
      
      # ======================================================
      # CONSERVATIVE PLOIDY CALL
      # ======================================================
      
      estimated.ploidy <- NA_real_
      call <- "ambiguous"
      
      
      # User-supplied known ploidy takes precedence
      if (
        !is.null(known.ploidy) &&
        current.group %in% names(known.ploidy)
      ) {
        
        estimated.ploidy <-
          known.ploidy[current.group]
        
        call <- "known"
        
      } else if (
        
        support.3n >= 0.20 &&
        support.3n >
        support.4n.outer + 0.10
        
      ) {
        
        estimated.ploidy <- 3
        call <- "triploid"
        
      } else if (
        
        support.4n.outer >= 0.15
        
      ) {
        
        estimated.ploidy <- 4
        call <- "tetraploid"
        
      } else if (
        
        support.2n >= 0.50 &&
        support.3n < 0.15 &&
        support.4n.outer < 0.10
        
      ) {
        
        # Important:
        # AF=0.5 can represent either AB diploid
        # or AABB balanced tetraploid.
        #
        # Do NOT automatically force CN2.
        estimated.ploidy <- NA_real_
        call <- "diploid_or_balanced_tetraploid"
        
      }
      
      
      data.frame(
        population = current.group,
        n.variants = n.variants,
        support.2n = support.2n,
        support.3n = support.3n,
        support.4n.outer = support.4n.outer,
        estimated.ploidy =
          estimated.ploidy,
        call = call,
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- dplyr::bind_rows(
    results
  )
  
  return(results)
}

.applyPloidyScaling <- function(
    TapestriExperiment,
    ploidy.table,
    sample.feature = "cluster",
    reference.ploidy = 2
) {
  
  groups <- as.character(
    SummarizedExperiment::colData(
      TapestriExperiment
    )[[sample.feature]]
  )
  
  lookup <- ploidy.table$estimated.ploidy
  
  names(lookup) <-
    ploidy.table$population
  
  
  cell.ploidy <- lookup[
    groups
  ]
  
  
  # populations for which ploidy could not be confidently
  # estimated remain unscaled
  cell.ploidy[
    is.na(cell.ploidy)
  ] <- reference.ploidy
  
  
  scaling.factor <-
    cell.ploidy /
    reference.ploidy
  
  
  SummarizedExperiment::colData(
    TapestriExperiment
  )$estimatedPloidy <- cell.ploidy
  
  
  SummarizedExperiment::colData(
    TapestriExperiment
  )$ploidyScalingFactor <-
    scaling.factor
  
  
  alt.exps <- c(
    "smoothedCopyNumberByChr",
    "smoothedCopyNumberByArm",
    "smoothedCopyNumberByCytob"
  )
  
  
  for (
    current.exp in
    intersect(
      alt.exps,
      SingleCellExperiment::altExpNames(
        TapestriExperiment
      )
    )
  ) {
    
    current.alt <-
      SingleCellExperiment::altExp(
        TapestriExperiment,
        current.exp
      )
    
    
    raw.cn <-
      SummarizedExperiment::assay(
        current.alt,
        "smoothedCopyNumber"
      )
    
    
    adjusted.cn <- sweep(
      raw.cn,
      MARGIN = 2,
      STATS = scaling.factor,
      FUN = "*"
    )
    
    
    SummarizedExperiment::assay(
      current.alt,
      "ploidyAdjustedSmoothedCopyNumber"
    ) <- adjusted.cn
    
    
    SummarizedExperiment::assay(
      current.alt,
      "ploidyAdjustedDiscreteCopyNumber"
    ) <- round(
      adjusted.cn
    )
    
    
    SingleCellExperiment::altExp(
      TapestriExperiment,
      current.exp
    ) <- current.alt
  }
  
  
  S4Vectors::metadata(
    TapestriExperiment
  )$ploidy.estimation <-
    ploidy.table
  
  
  return(
    TapestriExperiment
  )
}
