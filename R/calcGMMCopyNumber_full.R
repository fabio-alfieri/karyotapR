calcGMMCopyNumber.full <- function(TapestriExperiment,
                              cell.barcodes,
                              control.copy.number,
                              model.components = 1:5,
                              model.priors = NULL,
                              ...) {
  if (is.null(model.priors)) {
    model.priors <- rep(1, length(model.components))
  } else {
    if (length(model.priors) != length(model.components)) {
      cli::cli_abort("model.priors must be same length as model.components, or `NULL` for equal priors.")
    }
  }

  if (rlang::is_missing(control.copy.number)) {
    cli::cli_abort("{.arg control.copy.number} has not been set. Use {.fun karyotapR::generateControlCopyNumberTemplate}.")
  }

  if (length(cell.barcodes) == 0) {
    cli::cli_abort("cell.barcodes is empty.")
  } else {
    cli::cli_alert_info("Calculating GMMs using {length(cell.barcodes)} input cells.")
    filtered.tapestri.exp <- TapestriExperiment[, cell.barcodes]
  }

  # simulate probe counts
  simulated.norm.counts <- .generateSimulatedCNVCells(
    TapestriExperiment = filtered.tapestri.exp,
    control.copy.number = control.copy.number,
    ...
  )

  # smooth counts from simulated cells into smoothed copy number values
  cli::cli_progress_step("Fitting Gaussian distributions to simulated cells...")
  smoothing.method <- S4Vectors::metadata(TapestriExperiment)$smoothing.method
  smoothing.weights <- S4Vectors::metadata(TapestriExperiment)$smoothing.weights
  simulated.tapestri.experiment <- .smoothSimulatedCells(
    normalized.counts = simulated.norm.counts,
    probe.metadata = rowData(TapestriExperiment),
    smoothing.method = smoothing.method,
    smoothing.weights = smoothing.weights,
    ...
  )

  # fit Gaussian distributions to simulated cells
  cn.model.params.chr <- .fitGaussianDistributions(simulated.tapestri.experiment = simulated.tapestri.experiment, chromosome.scope = "chr")
  cn.model.params.arm <- .fitGaussianDistributions(simulated.tapestri.experiment = simulated.tapestri.experiment, chromosome.scope = "arm")
  cn.model.params.cytob <- .fitGaussianDistributions(simulated.tapestri.experiment = simulated.tapestri.experiment, chromosome.scope = "cytoband")

  # calculate posterior probabilities for each data point under each model component
  cli::cli_progress_step("Calculating posterior probabilities...")
  cn.model.table.chr <- .calcClassPosteriors(
    TapestriExperiment = TapestriExperiment,
    cn.model.params = cn.model.params.chr,
    model.components = model.components,
    model.priors = model.priors,
    chromosome.scope = "chr"
  )
  cn.model.table.arm <- .calcClassPosteriors(
    TapestriExperiment = TapestriExperiment,
    cn.model.params = cn.model.params.arm,
    model.components = model.components,
    model.priors = model.priors,
    chromosome.scope = "arm"
  )
  cn.model.table.cytob <- .calcClassPosteriors(
    TapestriExperiment = TapestriExperiment,
    cn.model.params = cn.model.params.cytob,
    model.components = model.components,
    model.priors = model.priors,
    chromosome.scope = "cytoband"
  )


  # call copy number values from posterior probabilities
  cli::cli_progress_step("Calling copy number from posterior probabilities...")
  cn.model.table.chr <- .callCopyNumberClasses(cn.model.table.chr)
  cn.model.table.arm <- .callCopyNumberClasses(cn.model.table.arm)
  cn.model.table.cytob <- .callCopyNumberClasses(cn.model.table.cytob)
  cli::cli_progress_done()

  # transform copy number calls to matrix
  # add copy number calls and model metadata to TapestriExperiment

  # whole chromosomes
  cli::cli_bullets(c("v" = "Saving whole chromosome copy number calls to altExp: smoothedCopyNumberByChr, assay: gmmCopyNumber..."))

  class.labels.chr.df <- cn.model.table.chr %>%
    dplyr::pull("cn.class") %>%
    purrr::map(\(x) tidyr::pivot_wider(x,
      names_from = "cell.barcode",
      values_from = "cn.class"
    )) %>%
    purrr::list_rbind() %>%
    as.data.frame() %>%
    magrittr::set_rownames(cn.model.table.chr$feature.id)

  SummarizedExperiment::assay(altExp(TapestriExperiment, "smoothedCopyNumberByChr"), "gmmCopyNumber") <- class.labels.chr.df

  # arms
  cli::cli_bullets(c("v" = "Saving chromosome arm copy number calls to altExp: smoothedCopyNumberByArm, assay: gmmCopyNumber..."))

  class.labels.arm.df <- cn.model.table.arm %>%
    dplyr::pull("cn.class") %>%
    purrr::map(\(x) tidyr::pivot_wider(x,
      names_from = "cell.barcode",
      values_from = "cn.class"
    )) %>%
    purrr::list_rbind() %>%
    as.data.frame() %>%
    magrittr::set_rownames(cn.model.table.arm$feature.id)

  SummarizedExperiment::assay(altExp(TapestriExperiment, "smoothedCopyNumberByArm"), "gmmCopyNumber") <- class.labels.arm.df

  # cytobands
  cli::cli_bullets(c("v" = "Saving chromosome arm copy number calls to altExp: smoothedCopyNumberByCytob, assay: gmmCopyNumber..."))

  class.labels.cytob.df <- cn.model.table.cytob %>%
    dplyr::pull("cn.class") %>%
    purrr::map(\(x) tidyr::pivot_wider(x,
      names_from = "cell.barcode",
      values_from = "cn.class"
    )) %>%
    purrr::list_rbind() %>%
    as.data.frame() %>%
    magrittr::set_rownames(cn.model.table.cytob$feature.id)

  SummarizedExperiment::assay(altExp(TapestriExperiment, "smoothedCopyNumberByCytob"), "gmmCopyNumber") <- class.labels.cytob.df

  TapestriExperiment@gmmParams <- list("chr" = cn.model.table.chr, "arm" = cn.model.table.arm, "cytoband" = cn.model.table.cytob)
  cli::cli_bullets(c("v" = "Saving GMM models and metadata to {.var gmmParams} slot..."))
  cli::cli_progress_done()

  return(TapestriExperiment)
}
