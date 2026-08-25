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
    min.variants = 50,
    min.variants.per.arm = 8,
    min.informative.arms = 10,
    tolerance = 0.05,
    min.triploid.arm.fraction = 0.55,
    known.ploidy = NULL
) {
  
  if (
    !"alleleFrequency" %in%
    SingleCellExperiment::altExpNames(TapestriExperiment)
  ) {
    cli::cli_abort(
      "alleleFrequency altExp not found."
    )
  }
  
  
  # ==========================================================
  # AF DATA
  # ==========================================================
  
  af.exp <- SingleCellExperiment::altExp(
    TapestriExperiment,
    "alleleFrequency"
  )
  
  af <- SummarizedExperiment::assay(
    af.exp,
    "alleleFrequency"
  )
  
  if (max(af, na.rm = TRUE) > 1.5) {
    af <- af / 100
  }
  
  
  variant.meta <- as.data.frame(
    SummarizedExperiment::rowData(
      af.exp
    )
  )
  
  
  # ==========================================================
  # GET ARM INFORMATION
  #
  # Match variant amplicon -> main experiment amplicon -> arm
  # ==========================================================
  
  if (
    "amplicon.id" %in% colnames(variant.meta)
  ) {
    
    idx <- match(
      variant.meta$amplicon.id,
      rownames(TapestriExperiment)
    )
    
    variant.meta$arm <-
      as.character(
        SummarizedExperiment::rowData(
          TapestriExperiment
        )$arm[idx]
      )
    
  } else {
    
    variant.meta$arm <- NA_character_
    
  }
  
  
  # ==========================================================
  # POPULATION LABELS
  # ==========================================================
  
  groups <- as.character(
    SummarizedExperiment::colData(
      TapestriExperiment
    )[[sample.feature]]
  )
  
  names(groups) <-
    colnames(TapestriExperiment)
  
  groups <- groups[
    colnames(af)
  ]
  
  
  group.names <- unique(groups)
  
  
  # ==========================================================
  # ANALYZE EACH POPULATION
  # ==========================================================
  
  results <- lapply(
    group.names,
    function(current.group) {
      
      
      # --------------------------------------------------------
      # User-defined ploidy takes precedence
      # --------------------------------------------------------
      
      if (
        !is.null(known.ploidy) &&
        current.group %in%
        names(known.ploidy)
      ) {
        
        return(
          tibble::tibble(
            population = current.group,
            n.variants = NA_integer_,
            n.informative.arms = NA_integer_,
            support.2n = NA_real_,
            support.3n = NA_real_,
            support.3n.left = NA_real_,
            support.3n.right = NA_real_,
            triploid.arm.fraction = NA_real_,
            estimated.ploidy =
              as.numeric(
                known.ploidy[
                  current.group
                ]
              ),
            call = "known"
          )
        )
      }
      
      
      cells <- which(
        groups == current.group
      )
      
      x <- af[
        ,
        cells,
        drop = FALSE
      ]
      
      
      # --------------------------------------------------------
      # Number of cells measured per SNP
      # --------------------------------------------------------
      
      n.obs <- rowSums(
        !is.na(x)
      )
      
      
      # --------------------------------------------------------
      # Pseudo-bulk VAF
      # --------------------------------------------------------
      
      pseudo.bulk.af <-
        matrixStats::rowMedians(
          x,
          na.rm = TRUE
        )
      
      
      df <- tibble::tibble(
        AF = pseudo.bulk.af,
        n.cells = n.obs,
        arm = variant.meta$arm
      ) %>%
        
        filter(
          n.cells >= min.cells,
          AF >= 0.10,
          AF <= 0.90,
          !is.na(arm)
        )
      
      
      n.variants <- nrow(df)
      
      
      # ======================================================
      # LOW INFORMATION
      # ======================================================
      
      if (
        n.variants < min.variants
      ) {
        
        return(
          tibble::tibble(
            population = current.group,
            n.variants = n.variants,
            n.informative.arms = NA_integer_,
            support.2n = NA_real_,
            support.3n = NA_real_,
            support.3n.left = NA_real_,
            support.3n.right = NA_real_,
            triploid.arm.fraction = NA_real_,
            estimated.ploidy = NA_real_,
            call = "low_information"
          )
        )
      }
      
      
      # ======================================================
      # GLOBAL SUPPORT
      # ======================================================
      
      support.2n <- mean(
        abs(
          df$AF - 0.5
        ) <= tolerance
      )
      
      
      support.3n.left <- mean(
        abs(
          df$AF - 1/3
        ) <= tolerance
      )
      
      
      support.3n.right <- mean(
        abs(
          df$AF - 2/3
        ) <= tolerance
      )
      
      
      support.3n <-
        support.3n.left +
        support.3n.right
      
      
      # ======================================================
      # ARM-LEVEL SUPPORT
      #
      # Each chromosome arm votes independently.
      # This avoids heavily represented arms dominating
      # the entire genome-wide ploidy call.
      # ======================================================
      
      arm.summary <- df %>%
        
        group_by(
          arm
        ) %>%
        
        filter(
          n() >=
            min.variants.per.arm
        ) %>%
        
        summarise(
          
          n.snps = n(),
          
          support2 = mean(
            abs(AF - 0.5)
            <= tolerance
          ),
          
          support3.left = mean(
            abs(AF - 1/3)
            <= tolerance
          ),
          
          support3.right = mean(
            abs(AF - 2/3)
            <= tolerance
          ),
          
          support3 =
            support3.left +
            support3.right,
          
          .groups = "drop"
          
        ) %>%
        
        mutate(
          
          state = case_when(
            
            support3 >
              support2 + 0.10 ~
              
              "3n-like",
            
            support2 >
              support3 + 0.10 ~
              
              "2n-like",
            
            TRUE ~
              
              "ambiguous"
          )
        )
      
      
      n.informative.arms <-
        nrow(
          arm.summary
        )
      
      
      if (
        n.informative.arms <
        min.informative.arms
      ) {
        
        return(
          tibble::tibble(
            population = current.group,
            n.variants = n.variants,
            n.informative.arms =
              n.informative.arms,
            support.2n = support.2n,
            support.3n = support.3n,
            support.3n.left =
              support.3n.left,
            support.3n.right =
              support.3n.right,
            triploid.arm.fraction =
              NA_real_,
            estimated.ploidy =
              NA_real_,
            call =
              "low_arm_information"
          )
        )
      }
      
      
      triploid.arm.fraction <- mean(
        arm.summary$state ==
          "3n-like"
      )
      
      
      diploid.arm.fraction <- mean(
        arm.summary$state ==
          "2n-like"
      )
      
      
      # ======================================================
      # STRICT TRIPLOID CALL
      # ======================================================
      
      triploid.call <-
        
        # enough information
        n.informative.arms >= min.informative.arms &&
        
        # reasonable genome-wide evidence
        support.3n >= 0.25 &&
        
        # both AAB and ABB signatures must exist
        support.3n.left >= 0.07 &&
        support.3n.right >= 0.07 &&
        
        # substantial fraction of arms must look triploid
        triploid.arm.fraction >= 0.45 &&
        
        # and triploid-like arms must outnumber diploid-like arms
        triploid.arm.fraction >=
        diploid.arm.fraction + 0.08
      
      # ======================================================
      # DIPLOID CALL
      # ======================================================
      
      diploid.call <-
        
        n.informative.arms >= min.informative.arms &&
        
        support.2n >= 0.30 &&
        
        diploid.arm.fraction >= 0.45 &&
        
        diploid.arm.fraction >=
        triploid.arm.fraction + 0.08
      
      
      # ======================================================
      # FINAL CLASSIFICATION
      # ======================================================
      
      if (triploid.call) {
        
        estimated.ploidy <- 3
        call <- "triploid"
        
      } else if (diploid.call) {
        
        estimated.ploidy <- 2
        call <- "diploid"
        
      } else {
        
        estimated.ploidy <- NA_real_
        call <- "ambiguous"
      }
      
      
      tibble::tibble(
        
        population =
          current.group,
        
        n.variants =
          n.variants,
        
        n.informative.arms =
          n.informative.arms,
        
        support.2n =
          support.2n,
        
        support.3n =
          support.3n,
        
        support.3n.left =
          support.3n.left,
        
        support.3n.right =
          support.3n.right,
        
        triploid.arm.fraction =
          triploid.arm.fraction,
        
        diploid.arm.fraction =
          diploid.arm.fraction,
        
        estimated.ploidy =
          estimated.ploidy,
        
        call =
          call
      )
    }
  )
  
  
  dplyr::bind_rows(
    results
  )
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
