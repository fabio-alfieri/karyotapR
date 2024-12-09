#' @param TapestriExperiment `TapestriExperiment` object.
#' @param copy.number Numeric, sets all entries of `copy.number` column in output. Default 2 (diploid).
#' @param sample.feature.label Character, sets all entries of `sample.label` column in output.
#'
#' @return `data.frame` with 3 columns named `arm`, `copy.number`, and `sample.label`
#' @export
#'
#' @describeIn calcCopyNumber generates a `data.frame` template for `control.copy.number` in `calcCopyNumber()`.
#' @order 2
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' control.copy.number <- generateControlCopyNumberTemplate(tap.object,
#'   copy.number = 2,
#'   sample.feature.label = "cellline1"
#' )
generateControlCopyNumberTemplate <- function(TapestriExperiment,
                                              copy.number = 2,
                                              sample.feature.label = NA) {
  if (any(is.na(unique(SummarizedExperiment::rowData(TapestriExperiment)$arm)))) {
    cli::cli_abort("Non-genomic probe found in rowData(<TapestriExperiment>)$arm column. Please remove before calculating copy number.")
  } 
  
  ploidy.template <- data.frame(
    cytoband = unique(paste0(SummarizedExperiment::rowData(TapestriExperiment)$chr,'_',
                             SummarizedExperiment::rowData(TapestriExperiment)$cytoband)),
    copy.number = copy.number,
    sample.label = sample.feature.label
  )
  rownames(ploidy.template) <- ploidy.template$cytoband
  
  return(ploidy.template)
}

#' @name calcCopyNumber
#'
#' @title Calculate relative copy number value for each cell-probe unit using reference sample
#'
#' @description `calcCopyNumber()` transforms the normalized count matrix `normcounts` of a `TapestriExperiment` object
#' into copy number values based on a set of reference cell barcodes and given copy number value (e.g. 2 for diploid).
#' This is practically used to set the median copy number of a usually diploid reference
#' cell population to a known copy number value, e.g. 2, and then calculate the copy number for all the
#' cells relative to that reference population. This occurs individually for each probe,
#' such that the result is one copy number value per cell barcode per probe (cell-probe unit).
#' `control.copy.number` is a `data.frame` lookup table used to indicate the copy number value and cell barcodes
#' to use as the reference. A template for `control.copy.number` can be generated using [generateControlCopyNumberTemplate()],
#' which will have a row for each chromosome arm represented in `TapestriExperiment`.
#'
#' The `control.copy.number` data.frame should include 3 columns named `arm`, `copy.number`, and `sample.label`.
#' `arm` is chromosome arm names from chr1p through chrXq, `copy.number` is the reference copy number value (2 = diploid), and `sample.label` is the
#' value corresponding to the `colData` column given in `sample.feature` to indicate the set of reference cell barcodes to use to set the copy number.
#' This is best used in a workflow where the cells are clustered first into their respective samples, and then one cluster is used as the reference population
#' the other clusters. This also allows for the baseline copy number to be set for each chromosome arm individually in the case where the
#' reference population is not completely diploid.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param control.copy.number `data.frame` with columns `arm`, `copy.number`, and `sample.label`. See details.
#' @param sample.feature Character, `colData` column to use for subsetting cell.barcodes. Default "cluster".
#' @param remove.bad.probes Logical, if `TRUE`, probes with median normalized counts = 0 are removed from the returned `TapestriExperiment`. If FALSE (default), probes with median normalized counts = 0 throw error and stop function.
#'
#' @return `TapestriExperiment` object with cell-probe copy number values in `copyNumber` assay slot.
#' @export
#'
#' @rdname calcCopyNumber
#' @order 1
#'
#' @concept copy number
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' tap.object <- calcNormCounts(tap.object)
#' control.copy.number <- generateControlCopyNumberTemplate(tap.object,
#'   copy.number = 2,
#'   sample.feature.label = "cellline1"
#' )
#' tap.object <- calcCopyNumber(tap.object,
#'   control.copy.number,
#'   sample.feature = "test.cluster"
#' )
calcCopyNumber <- function(TapestriExperiment,
                           control.copy.number,
                           sample.feature = "cluster",
                           remove.bad.probes = FALSE) {
  sample.feature <- tolower(sample.feature)
  
  # error checks
  if (!sample.feature %in% colnames(SummarizedExperiment::colData(TapestriExperiment))) {
    cli::cli_abort("{.var sample.feature} {.q {sample.feature}} not found in {.var colData}.")
  }
  
  if (any(!unique(control.copy.number$sample.label) %in% unique(SummarizedExperiment::colData(TapestriExperiment)[, sample.feature]))) {
    cli::cli_abort("{.var control.copy.number} {.q sample.label} elements not found in {colData.} Check {.var control.copy.number.}")
  }
  
  counts.mat <- SummarizedExperiment::assay(TapestriExperiment, "normcounts")
  
  # get median normalized counts for each probe based on control.copy.number
  probe.table <- as.data.frame(SummarizedExperiment::rowData(TapestriExperiment))[, c("probe.id", "chr", "cytoband")]
  probe.table$cytoband <- paste0(probe.table$chr,probe.table$cytoband)
  probe.table$chr <- NULL
  probe.table <- merge(probe.table, control.copy.number, by = "cytoband", all.x = TRUE, sort = FALSE)
  rownames(probe.table) <- probe.table$probe.id
  
  sample.feature.lookup <- SummarizedExperiment::colData(TapestriExperiment)[, sample.feature, drop = FALSE]
  
  # define function for calculating median from cell subset
  getProbeMedian <- function(idx) {
    probe.info <- probe.table[idx, ]
    probe.median <- median(counts.mat[probe.info$probe.id, rownames(sample.feature.lookup)[sample.feature.lookup[, 1] == probe.info$sample.label]])
    return(probe.median)
  }
  
  probe.medians <- lapply(seq_len(nrow(probe.table)), getProbeMedian)
  probe.medians <- unlist(probe.medians)
  names(probe.medians) <- probe.table$probe.id
  
  
  # check for probes with median = 0
  bad.probes <- NULL
  if (any(probe.medians == 0)) {
    if (remove.bad.probes == FALSE) {
      cli::cli_abort("{names(probe.medians[probe.medians == 0])} control cell median equal to 0. Filter out prior to proceeding.")
    } else {
      bad.probes <- names(probe.medians)[which(probe.medians == 0)]
    }
  }
  
  probe.medians <- probe.medians[rownames(SummarizedExperiment::rowData(TapestriExperiment))] # reorder based on rowData
  counts.ploidy <- sweep(x = counts.mat, 1, probe.medians, "/") # normalize relative to medians
  probe.table <- probe.table[rownames(SummarizedExperiment::rowData(TapestriExperiment)), ] # reorder based on rowData
  counts.ploidy <- sweep(x = counts.ploidy, 1, probe.table$copy.number, "*") # scale to control copy number
  
  SummarizedExperiment::assay(TapestriExperiment, "copyNumber") <- counts.ploidy
  
  if (!is.null(bad.probes)) {
    TapestriExperiment <- TapestriExperiment[setdiff(rownames(TapestriExperiment), bad.probes), ]
    cli::cli_alert_info("Probes removed for 0 median value: {.q {bad.probes}}.")
  }
  
  return(TapestriExperiment)
}



#' Smooth copy number values across chromosomes and chromosome arms
#'
#' `calcSmoothCopyNumber()` takes `copyNumber` slot values for probes on a chromosome and smooths them by median (default) for each chromosome
#' and chromosome arm, resulting in one copy number value per chromosome and chromosome arm for each cell barcode.
#' Cell-chromosome values are then discretized into integers by conventional rounding (1.5 <= x < 2.5 rounds to 2).
#' Smoothed copy number and discretized smoothed copy number values are stored as `smoothedCopyNumber` and `discreteCopyNumber` assays,
#' in `altExp` slots `smoothedCopyNumberByChr` for chromosome-level smoothing, and `smoothedCopyNumberByArm` for chromosome arm-level smoothing.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param method Character, smoothing method: median (default) or mean.
#'
#' @importFrom rlang .data
#'
#' @return `TapestriExperiment` with `smoothedCopyNumber` and `discreteCopyNumber` assays in `altExp` slots `smoothedCopyNumberByChr` and `smoothedCopyNumberByArm`.
#' @export
#'
#' @concept copy number
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' tap.object <- calcNormCounts(tap.object)
#' control.copy.number <- generateControlCopyNumberTemplate(tap.object,
#'   copy.number = 2,
#'   sample.feature.label = "cellline1"
#' )
#' tap.object <- calcCopyNumber(tap.object,
#'   control.copy.number,
#'   sample.feature = "test.cluster"
#' ) 
#' tap.object <- calcSmoothCopyNumber(tap.object) 
calcSmoothCopyNumber <- function(TapestriExperiment, 
                                 method = "median") {
  method <- tolower(method)
  
  if (method == "median") {
    smooth.func <- stats::median
  } else if (method == "mean") {
    smooth.func <- mean
  } else {
    cli::cli_abort("{.var method} {.q {method}}, not recognized. Please use {.q mean} or {.q median}.")
  }
  
  if(!"copyNumber" %in% SummarizedExperiment::assayNames(TapestriExperiment)){
    cli::cli_abort("{.q copyNumber} assay not found in {.code TapestriExperiment} object. Did you run {.fn karyotapR::calcCopyNumber} first?")
  }
  
  cli::cli_progress_step("Smoothing copy number by {method}...", )
  
  ploidy.counts <- SummarizedExperiment::assay(TapestriExperiment, "copyNumber")
  
  ploidy.tidy <- ploidy.counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("probe.id") %>%
    tidyr::pivot_longer(
      cols = !tidyr::matches("probe.id"),
      names_to = "cell.barcode",
      values_to = "ploidy"
    ) %>%
    dplyr::left_join(as.data.frame(SummarizedExperiment::rowData(TapestriExperiment)[, c("probe.id", "chr", "arm", "cytoband")]), by = "probe.id")
  
  smoothed.ploidy.chr <- ploidy.tidy %>%
    dplyr::group_by(.data$cell.barcode, .data$chr) %>%
    dplyr::summarize(
      smooth.ploidy = smooth.func(.data$ploidy),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of("chr"),
      values_from = dplyr::all_of("smooth.ploidy"),
      names_from = dplyr::all_of("cell.barcode")
    ) %>%
    tibble::column_to_rownames("chr")
  
  # reorder to match input matrix
  smoothed.ploidy.chr <- smoothed.ploidy.chr[, colnames(ploidy.counts)]
  
  smoothed.ploidy.arm <- ploidy.tidy %>%
    dplyr::group_by(.data$cell.barcode, .data$arm) %>%
    dplyr::summarize(
      smooth.ploidy = smooth.func(.data$ploidy),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of("arm"),
      values_from = dplyr::all_of("smooth.ploidy"),
      names_from = dplyr::all_of("cell.barcode")
    ) %>%
    tibble::column_to_rownames("arm")
  
  # reorder to match input matrix
  smoothed.ploidy.arm <- smoothed.ploidy.arm[, colnames(ploidy.counts)]
  
  smoothed.ploidy.cytob <- ploidy.tidy %>%
    dplyr::group_by(.data$cell.barcode, .data$cytoband) %>%
    dplyr::summarize(
      smooth.ploidy = smooth.func(.data$ploidy),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of("cytoband"),
      values_from = dplyr::all_of("smooth.ploidy"),
      names_from = dplyr::all_of("cell.barcode")
    ) %>%
    tibble::column_to_rownames("cytoband")
  
  discrete.ploidy.chr <- round(smoothed.ploidy.chr, 0)
  discrete.ploidy.arm <- round(smoothed.ploidy.arm, 0)
  discrete.ploidy.cytob <- round(smoothed.ploidy.cytob, 0)
  
  
  smoothed.ploidy.chr <- SingleCellExperiment::SingleCellExperiment(list(
    smoothedCopyNumber = smoothed.ploidy.chr,
    discreteCopyNumber = discrete.ploidy.chr
  ))
  
  smoothed.ploidy.arm <- SingleCellExperiment::SingleCellExperiment(list(
    smoothedCopyNumber = smoothed.ploidy.arm,
    discreteCopyNumber = discrete.ploidy.arm
  ))
  
  smoothed.ploidy.cytob <- SingleCellExperiment::SingleCellExperiment(list(
    smoothedCopyNumber = smoothed.ploidy.cytob,
    discreteCopyNumber = discrete.ploidy.cytob
  ))
  
  smoothed.ploidy.chr <- .TapestriExperiment(smoothed.ploidy.chr)
  smoothed.ploidy.arm <- .TapestriExperiment(smoothed.ploidy.arm)
  smoothed.ploidy.cytob <- .TapestriExperiment(smoothed.ploidy.cytob)
  
  SingleCellExperiment::altExp(TapestriExperiment, "smoothedCopyNumberByChr", withDimnames = TRUE) <- smoothed.ploidy.chr
  SingleCellExperiment::altExp(TapestriExperiment, "smoothedCopyNumberByArm", withDimnames = TRUE) <- smoothed.ploidy.arm
  SingleCellExperiment::altExp(TapestriExperiment, "smoothedCopyNumberByCytob", withDimnames = TRUE) <- smoothed.ploidy.cytob
  
  return(TapestriExperiment)
}