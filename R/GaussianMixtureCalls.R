#' Call copy number for each cell-chromosome using Gaussian mixture models
#'
#' Uses control cells to simulate expected smoothed copy number distributions for all chromosomes across each of `model.components` (copy number level).
#' Then uses the distributions to calculate posterior probabilities for each cell-chromosome belonging to each of copy number level.
#' Each cell-chromosome is assigned the copy number value for which its posterior probability is highest.
#' This is done for both whole chromosomes and chromosome arms.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param cell.barcodes character, vector of cell barcodes to fit GMM. Usually corresponds to diploid control.
#' @param control.copy.number `data.frame` with columns `arm` and `copy.number`, indicating of known copy number of cells in `cell.barcodes`.
#' @param model.components numeric, vector of copy number GMM components to calculate, default `1:5` (for copy number = 1, 2, 3, 4, 5).
#' @param model.priors numeric, relative prior probabilities for each GMM component. If `NULL` (default), assumes equal priors.
#' @param ... Additional parameters to be passed to internal functions.
#'
#' @return `TapestriExperiment` object with copy number calls based on the calculated GMMs, saved to `gmmCopyNumber` slot of `smoothedCopyNumberByChr` and `smoothedCopyNumberByArm` altExps.
#' GMM parameters for each `feature.id` are saved to the `metadata` slot.
#' @export
#'
#' @concept copy number
#'
#' @examples
#' \dontrun{
#' control.copy.number <- generateControlCopyNumberTemplate(TapestriExperiment, 2)
#' TapestriExperiment <- calcGMMCopyNumber(TapestriExperiment,
#'   cell.barcodes = colnames(TapestriExperiment),
#'   control.copy.number = control.copy.number,
#'   model.components = 1:5,
#'   model.priors = c(1, 1, 1, 1, 1)
#' )
#' }
calcGMMCopyNumber <- function(TapestriExperiment,
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
    cn.model.params = cn.model.params.arm,
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
  cli::cli_bullets(c("v" = "Saving chromosome arm copy number calls to altExp: smoothedCopyNumberByArm, assay: gmmCopyNumber..."))

  class.labels.cytob.df <- cn.model.tabl.cytob %>%
    dplyr::pull("cn.class") %>%
    purrr::map(\(x) tidyr::pivot_wider(x,
      names_from = "cell.barcode",
      values_from = "cn.class"
    )) %>%
    purrr::list_rbind() %>%
    as.data.frame() %>%
    magrittr::set_rownames(cn.model.tabl.cytob$feature.id)

  SummarizedExperiment::assay(altExp(TapestriExperiment, "smoothedCopyNumberByArm"), "gmmCopyNumber") <- class.labels.cytob.df

  TapestriExperiment@gmmParams <- list("chr" = cn.model.table.chr, "arm" = cn.model.table.arm, "cytoband" = cn.model.table.cytob)
  cli::cli_bullets(c("v" = "Saving GMM models and metadata to {.var gmmParams} slot..."))
  cli::cli_progress_done()

  return(TapestriExperiment)
}

# model probes and generate simulated cells with normalized counts
.generateSimulatedCNVCells <- function(TapestriExperiment,
                                       control.copy.number,
                                       n.simulated.cells = 500) {
  cli::cli_progress_step("Generating probe values for {n.simulated.cells} simulated cells...")

  raw.counts <- SummarizedExperiment::assay(TapestriExperiment, "counts")
  norm.counts <- .MBNormCounts(raw.counts)
  norm.counts[norm.counts == 0] <- 1 # pseudocount zeros to ones
  norm.counts <- as.list(as.data.frame(t(norm.counts))) # convert to list of probes

  # fit probe parameters
  probe.model.fit <- suppressWarnings(norm.counts %>% purrr::map(\(x) fitdistrplus::fitdist(data = x, distr = "weibull")))
  probe.model.fit <- probe.model.fit %>%
    purrr::map(~ .x$estimate) %>%
    purrr::list_transpose() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("probe.id")

  # get control copy number, calculate scalar values for copy number = 1:6,
  # multiply scale parameter by scalar value
  probe.transform <- SingleCellExperiment::rowData(TapestriExperiment) %>%
    tibble::as_tibble() %>%
    dplyr::select("probe.id","chr", "cytoband") 
  probe.transform$cytoband <- paste0('chr',probe.transform$chr,probe.transform$cytoband)
  probe.transform <- probe.transform %>%
    dplyr::inner_join(control.copy.number %>% dplyr::select(!c("sample.label")), by = "cytoband")

  probe.transform <- probe.transform %>% dplyr::mutate("scalar" = purrr::map(.data$copy.number, \(x) 1:6 / x))
  probe.transform <- dplyr::inner_join(probe.transform, probe.model.fit, by = "probe.id")
  probe.transform <- probe.transform %>% dplyr::mutate("scale.transform" = purrr::map2(.data$scalar, .data$scale, \(x, y) x * y))
  probe.transform$probe.id <- factor(probe.transform$probe.id, levels = unique(probe.transform$probe.id))

  # simulate normalized counts for each parameter set
  simulated.counts <- probe.transform %>%
    tidyr::unnest("scale.transform") %>%
    dplyr::mutate("sim.counts" = purrr::map2(
      .data$shape,
      .data$scale.transform,
      function(shape, scale.transform) {
        stats::rweibull(
          n = n.simulated.cells,
          shape = shape, scale =
            scale.transform
        )
      }
    ))

  # combine counts into matrix
  simulated.counts <- simulated.counts %>% dplyr::select("probe.id", "sim.counts")
  simulated.counts <- split(simulated.counts, simulated.counts$probe.id) %>%
    purrr::map(\(x) unlist(x$sim.counts)) %>%
    as.data.frame() %>%
    t()
  colnames(simulated.counts) <- paste0("sim_cn", 1:6) %>%
    purrr::map(\(x) paste(x, seq_len(n.simulated.cells), sep = "_")) %>%
    unlist()

  return(simulated.counts)
}

# smooth simulated counts into copy number values
.smoothSimulatedCells <- function(normalized.counts,
                                  probe.metadata,
                                  smoothing.method,
                                  smoothing.weights,
                                  genome = "hg19") {
  # generate tapestri experiment object and copy normcounts slot
  tapestri.sim <- .createTapestriExperiment.sim(
    counts = normalized.counts,
    probe.metadata = probe.metadata,
    genome = genome
  )
  SummarizedExperiment::assay(tapestri.sim, "normcounts") <- normalized.counts

  # delimit cell barcodes to get copy number classes
  cn.sim.class <- matrix(unlist(strsplit(colData(tapestri.sim)$cell.barcode, split = "_")), ncol = 3, byrow = TRUE)
  cn.sim.class <- paste(cn.sim.class[, 1], cn.sim.class[, 2], sep = "_")

  SummarizedExperiment::colData(tapestri.sim)$cn.sim.class <- as.factor(cn.sim.class)

  control.cn <- generateControlCopyNumberTemplate(tapestri.sim, copy.number = 2, sample.feature.label = "sim_cn2")
  control.cn$cytoband <- str_replace(control.cn$cytoband, '_','')
  tapestri.sim <- calcCopyNumber(tapestri.sim, control.copy.number = control.cn, sample.feature = "cn.sim.class")
  
  cli::cli_progress_step("Smoothing copy number by {smoothing.method}...")
  tapestri.sim <- suppressMessages(calcSmoothCopyNumber(tapestri.sim, method = smoothing.method, 
                                   control.copy.number = control.cn, sample.feature = "cn.sim.class", weight.range = smoothing.weights))
}

# fit Gaussian distributions to simulated cells
.fitGaussianDistributions <- function(simulated.tapestri.experiment,
                                      chromosome.scope) {
  if (chromosome.scope == "chr" | chromosome.scope == "chromosome") {
    sim.data.tidy <- getTidyData(simulated.tapestri.experiment,
      alt.exp = "smoothedCopyNumberByChr",
      assay = "smoothedCopyNumber"
    )
  } else if (chromosome.scope == "arm") {
    sim.data.tidy <- getTidyData(simulated.tapestri.experiment,
      alt.exp = "smoothedCopyNumberByArm",
      assay = "smoothedCopyNumber"
    )
  } else if (chromosome.scope == "cytoband") {
    sim.data.tidy <- getTidyData(simulated.tapestri.experiment,
      alt.exp = "smoothedCopyNumberByCytob",
      assay = "smoothedCopyNumber"
    )
  } else {
    cli::cli_abort("chromosome.scope should be 'chr', 'arm' or 'cytoband")
  }

  # fit Gaussian distributions
  cn.model.params <- sim.data.tidy %>%
    tidyr::nest(.by = c("feature.id", "cn.sim.class"), .key = "smoothed.counts") %>%
    dplyr::mutate("fit" = purrr::map(.data$smoothed.counts, \(x) suppressWarnings(fitdistrplus::fitdist(x$smoothedCopyNumber, distr = "norm")), .progress = TRUE)) %>%
    dplyr::mutate("params" = purrr::map(.data$fit, \(x) x$estimate)) %>%
    tidyr::unnest_wider("params") %>%
    dplyr::select("feature.id", "cn.sim.class", "mean", "sd")

  return(cn.model.params)
}

# get probabilities for belonging to a given component of the Gaussian mixture
.calcClassPosteriors <- function(TapestriExperiment, cn.model.params, model.components, model.priors, chromosome.scope) {
  components.filtered <- paste0("sim_cn", model.components)

  if (chromosome.scope == "chr" | chromosome.scope == "chromosome") {
    sim.data.tidy <- getTidyData(TapestriExperiment,
      alt.exp = "smoothedCopyNumberByChr",
      assay = "smoothedCopyNumber"
    )
  } else if (chromosome.scope == "arm") {
    sim.data.tidy <- getTidyData(TapestriExperiment,
      alt.exp = "smoothedCopyNumberByArm",
      assay = "smoothedCopyNumber"
    )
  } else if (chromosome.scope == "cytoband") {
    sim.data.tidy <- getTidyData(TapestriExperiment,
      alt.exp = "smoothedCopyNumberByCytob",
      assay = "smoothedCopyNumber"
    )
  } else {
    cli::cli_abort("chromosome.scope should be 'chr', 'arm' or 'cytoband'")
  }

  # get smoothed copy number and combine in tibble with copy number model parameters
  smoothed.cn.df <- sim.data.tidy %>%
    dplyr::select("feature.id", "cell.barcode", "smoothedCopyNumber") %>%
    tidyr::nest(.by = "feature.id", .key = "smoothed.cn")
  cn.params.df <- cn.model.params %>%
    dplyr::filter(.data$cn.sim.class %in% components.filtered) %>%
    tidyr::nest(.by = "feature.id", .key = "model")
  gmm.table <- dplyr::inner_join(smoothed.cn.df, cn.params.df, by = "feature.id")

  # iterate over both lists in parallel to get probability density function values across all GMM components for all data points
  gmm.table <- gmm.table %>% dplyr::mutate("pdf" = purrr::map2(.data$smoothed.cn, .data$model, function(smoothed.cn, model) {
    df <- purrr::pmap(model, function(cn.sim.class, mean, sd) {
      stats::setNames(data.frame(stats::dnorm(smoothed.cn[["smoothedCopyNumber"]], mean = mean, sd = sd)), cn.sim.class)
    }) %>%
      purrr::list_cbind()
    rownames(df) <- smoothed.cn[["cell.barcode"]]
    df <- tibble::as_tibble(df)
    return(df)
  }))

  # get Bayes theorem denominator for each GMM (aka evidence or marginal likelihood)
  gmm.table <- gmm.table %>%
    dplyr::mutate("model.evidence" = purrr::map(.data$pdf, function(pdf) {
      as.matrix(pdf) %*% model.priors %>% drop()
    }))

  # calculate probability of a data point belonging to each copy number class (i.e. copy number class given a data point)
  gmm.table <- gmm.table %>% dplyr::mutate("cn.probability" = purrr::map2(.data$pdf, .data$model.evidence, function(pdf, model.evidence) {
    probs <- sweep(pdf, 2, model.priors, "*") # priors multiplied through each column
    probs <- sweep(probs, 1, model.evidence, "/") # model.evidence divided through each row
    probs <- tibble::as_tibble(probs)
    return(probs)
  }))

  return(gmm.table)
}

# assign copy numbers to data points by highest posterior probability
.callCopyNumberClasses <- function(cn.model.table) {
  # iterate over cn.probability and smoothed.cn lists in parallel
  cn.model.table <- cn.model.table %>%
    dplyr::mutate(cn.class = purrr::map2(.data$cn.probability, .data$smoothed.cn, function(cn.probability, smoothed.cn) {
      result <- max.col(cn.probability, "first") # pick class with highest probability
      result <- colnames(cn.probability)[result] # label selections
      result <- as.numeric(matrix(unlist(strsplit(result, split = "cn")), ncol = 2, byrow = TRUE)[, 2]) # pull copy number value from class label
      result <- tibble::as_tibble(data.frame(cell.barcode = purrr::pluck(smoothed.cn, "cell.barcode"), cn.class = result))
    }))

  return(cn.model.table)
}

#' Calculate decision boundaries between components of copy number GMMs
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param chromosome.scope "chr" or "arm", for using models for either whole chromosomes or chromosome arms. Default "chr".
#'
#' @return tibble containing boundary values of GMMs for each `feature.id`.
#' @export
#'
#' @concept copy number
#'
#' @examples
#' \dontrun{
#' control.copy.number <- generateControlCopyNumberTemplate(TapestriExperiment, 2)
#' boundaries <- getGMMBoundaries(TapestriExperiment,
#'   chromosome.scope = "chr"
#' )
#' }
getGMMBoundaries <- function(TapestriExperiment, chromosome.scope = "chr") {
  if (S4Vectors::isEmpty(TapestriExperiment@gmmParams)) {
    cli::cli_abort("GMM metadata not found. Did you run `calcGMMCopyNumber()` yet?")
  }

  if (chromosome.scope == "chr" | chromosome.scope == "chromosome") {
    model.metadata <- TapestriExperiment@gmmParams$chr
  } else if (chromosome.scope == "arm") {
    model.metadata <- TapestriExperiment@gmmParams$arm
  } else if (chromosome.scope == "cytoband") {
    model.metadata <- TapestriExperiment@gmmParams$cytoband
  } else {
    cli::cli_abort("chromosome.scope should be 'chr or 'arm'")
  }

  model.boundaries <- model.metadata %>%
    dplyr::select("feature.id", "model") %>%
    dplyr::group_by(.data$feature.id) %>%
    dplyr::transmute("boundary" = purrr::map(.data$model, ~ .singleModelBoundary(.x))) %>%
    tidyr::unnest_wider(col = "boundary", names_sep = ".") %>%
    dplyr::ungroup()

  return(model.boundaries)
}

.getSingleBoundary <- function(mean1, sd1, mean2, sd2) {
  normals <- data.frame(
    idx = seq(0, 10, 0.01),
    component1 = stats::dnorm(x = seq(0, 10, 0.01), mean = mean1, sd = sd1),
    component2 = stats::dnorm(x = seq(0, 10, 0.01), mean = mean2, sd = sd2)
  )

  # calculate absolute minimum between two curves
  normals <- normals %>% dplyr::mutate(delta = abs(.data$component1 - .data$component2))

  # search for absolute minimum between the two means
  normals.filtered <- normals %>% dplyr::filter(.data$idx >= mean1, .data$idx <= mean2)
  boundary <- normals.filtered$idx[which(normals.filtered$delta == min(normals.filtered$delta))]

  return(boundary)
}

.singleModelBoundary <- function(current.model) {
  current.model.boundaries <- list()
  for (i in seq_len(nrow(current.model) - 1)) {
    current.model.boundaries[i] <- .getSingleBoundary(
      mean1 = current.model$mean[i], sd1 = current.model$sd[i],
      mean2 = current.model$mean[i + 1], sd2 = current.model$sd[i + 1]
    )
  }
  current.model.boundaries <- unlist(current.model.boundaries)
  return(current.model.boundaries)
}

#' Plot copy number GMM components
#'
#' Plots the probability densities of GMM components for given chromosome or chromosome arm, store in a `TapestriExperiment`.
#' [calcGMMCopyNumber()] must be run first.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param feature.id chromosome or chromosome arm to plot.
#' @param draw.boundaries logical, if `TRUE`, draw decision boundaries between each Gaussian component.
#' @param chromosome.scope "chr" or "arm", for plotting models for either whole chromosomes or chromosome arms.
#'
#' @return `ggplot` object
#' @export
#'
#' @import ggplot2
#'
#' @concept copy number
#' @concept plots
#'
#' @examples
#' \dontrun{
#' TapestriExperiment <- plotCopyNumberGMM(TapestriExperiment,
#'   feature.id = 7,
#'   chromosome.scope = "chr",
#'   draw.boundaries = TRUE
#' )
#' }
plotCopyNumberGMM <- function(TapestriExperiment,
                              feature.id = 1,
                              chromosome.scope = "chr",
                              draw.boundaries = FALSE) {
  if (S4Vectors::isEmpty(TapestriExperiment@gmmParams)) {
    cli::cli_abort("GMM metadata not found. Did you run `calcGMMCopyNumber()` yet?")
  }

  if (chromosome.scope == "chr" | chromosome.scope == "chromosome") {
    model.metadata <- TapestriExperiment@gmmParams$chr
  } else if (chromosome.scope == "arm") {
    model.metadata <- TapestriExperiment@gmmParams$arm
  } else if (chromosome.scope == "cytoband") {
    model.metadata <- TapestriExperiment@gmmParams$cytoband
  } else {
    cli::cli_abort("chromosome.scope should be 'chr', 'arm' or 'cytoband'")
  }

  if (length(feature.id) != 1) {
    cli::cli_abort("`feature.id` must have length = 1")
  }

  if (!feature.id %in% model.metadata$feature.id) {
    cli::cli_abort("`feature.id` not found in model metadata. Check `feature.id` and `chromosome.scope`.")
  }

  model.fit <- model.metadata[model.metadata[, "feature.id"] == feature.id, ] %>%
    purrr::pluck("model", 1) %>%
    purrr::pmap(function(mean, sd, ...) data.frame("y" = stats::dnorm(seq(0, 10, 0.05), "mean" = mean, "sd" = sd))) %>%
    purrr::map(\(x) tibble::add_column(x, "x" = seq(0, 10, 0.05), .before = 0)) %>%
    dplyr::bind_rows(.id = "cn")

  model.fit$y <- model.fit$y / max(model.fit$y) # normalize to max = 1

  model.plot <- model.fit %>%
    ggplot(aes(x = .data$x, y = .data$y, color = .data$cn)) +
    geom_line() +
    theme_bw() +
    labs(title = paste("Chromosome", feature.id), y = "Density", x = "Copy Number") +
    scale_x_continuous(breaks = 0:10)

  if (draw.boundaries == TRUE) {
    model.boundaries <- .singleModelBoundary(current.model = model.metadata[model.metadata[, "feature.id"] == feature.id, ] %>% purrr::pluck("model", 1))
    model.plot <- model.plot + geom_vline(xintercept = model.boundaries, linetype = "dashed")
  }

  return(model.plot)
}
