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
  
  ploidy.tidy$cytoband <- paste0(substr(ploidy.tidy$arm, 1, nchar(as.character(ploidy.tidy$arm))-1), ploidy.tidy$cytoband)

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
  
  # reorder to march input matrix
  smoothed.ploidy.cytob <- smoothed.ploidy.cytob[, colnames(ploidy.counts)]
  
  # reorder cytobands
  smoothed.ploidy.cytob <- smoothed.ploidy.cytob[order(match(rownames(smoothed.ploidy.cytob), 
               ploidy.tidy$cytoband)), , drop = FALSE]


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
