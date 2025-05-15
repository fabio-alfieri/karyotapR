#' Cluster assay data by Principal Components Analysis
#'
#' Analyzes assay data by Principal Components Analysis (PCA) and saves results
#' to `reducedDims` slot of `TapestriObject`.
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param alt.exp Character, `altExp` to use, `NULL` uses top-level/main experiment. Default "alleleFrequency".
#' @param assay Character, `assay` to use. `NULL` (default) uses first-indexed assay.
#' @param sd.min.threshold Numeric, minimum threshold for allelefreq.sd. Increase to run PCA on fewer, more variable dimensions. Set to `NULL` if not using for alleleFrequency slot. Default `NULL`.
#' @param center Logical, if `TRUE` (default), variables are shifted to be zero centered. See [stats::prcomp()].
#' @param scale. Logical,if `TRUE` (default), variables are scaled to have unit variance prior to PCA. See [stats::prcomp()].
#'
#' @return `TapestriExperiment` with PCA results saved to `reducedDims` slot of `altExp`, and proportion of variance explained by each PC saved to `metadata` slot of `altExp.`
#' @export
#'
#' @seealso [stats::prcomp()] for PCA method details.
#'
#' @concept dimensional reduction
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment
#' tap.object <- runPCA(tap.object, alt.exp = "alleleFrequency")
runPCA <- function(TapestriExperiment,
                   alt.exp = "alleleFrequency",
                   assay = NULL,
                   sd.min.threshold = NULL,
                   center = TRUE,
                   scale. = TRUE) {
  assay <- .SelectAssay(TapestriExperiment, alt.exp = alt.exp, assay = assay)

  # filter by sd.min.threshold
  if (!is.null(sd.min.threshold)) {
    feature.set.filter <- SingleCellExperiment::rowData(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$allelefreq.sd >= sd.min.threshold
    if (all(!feature.set.filter)) {
      cli::cli_abort("All features filtered out. Reduce sd.min.threshold.")
    } else {
      if (alt.exp == assay) {
        cli::cli_alert_info("Running PCA on {.emph {alt.exp}} with {sum(feature.set.filter)} features.")
      } else {
        cli::cli_alert_info("Running PCA on {.emph {alt.exp}} {.emph {assay}} with {sum(feature.set.filter)} features.")
      }
    }
  }
 
  pca.assay <- getTidyData(TapestriExperiment,
    alt.exp = alt.exp,
    assay = assay,
    feature.id.as.factor = FALSE
  )

  pca.assay <- pca.assay %>%
    tidyr::pivot_wider(
      id_cols = "feature.id",
      names_from = "cell.barcode",
      values_from = {{ assay }}
    ) %>%
    tibble::column_to_rownames("feature.id") %>%
    as.data.frame()

  if (!is.null(sd.min.threshold)) {
    pca.assay <- pca.assay[feature.set.filter, ]
  }

  # transpose matrix
  pca.assay <- t(pca.assay)

  # run PCA
  pca.result <- stats::prcomp(pca.assay, center = center, scale. = scale.)

  # save PCA to main experiment
  if (is.null(alt.exp)) {
    SingleCellExperiment::reducedDim(TapestriExperiment, "PCA", withDimnames = TRUE) <- pca.result$x
    S4Vectors::metadata(TapestriExperiment)$pca.proportion.of.variance <-
      summary(pca.result)$importance["Proportion of Variance", ]
    S4Vectors::metadata(TapestriExperiment)$pca.assay <- assay
  } else { # save PCA to alt experiment
    SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, alt.exp), "PCA", withDimnames = TRUE) <- pca.result$x
    S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.proportion.of.variance <-
      summary(pca.result)$importance["Proportion of Variance", ]
    S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.assay <- assay
  }

  return(TapestriExperiment)
}

#' Plot of PCA proportion of variance explained
#'
#' Draws "knee plot" of PCA proportion of variance explained to determine which principal components (PCs) to include for downstream applications e.g. clustering.
#' Variance explained for each PC is indicated by the line.
#' Cumulative variance explained is indicated by the bars.
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param alt.exp Character, `altExp` to use, `NULL` uses top-level/main experiment. Default "alleleFrequency".
#' @param n.pcs Numeric, number of PCs to plot, starting at 1. Default 10.
#'
#' @return ggplot2 object, combined line plot and bar graph
#' @export
#'
#' @concept dimensional reduction
#' @concept plots
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' tap.object <- runPCA(tap.object, alt.exp = "alleleFrequency")
#' PCAKneePlot(tap.object, n.pcs = 5)
PCAKneePlot <- function(TapestriExperiment, alt.exp = "alleleFrequency", n.pcs = 10) {
  if (is.null(alt.exp)) {
    knee.df <- data.frame("prop.variance" = S4Vectors::metadata(TapestriExperiment)$pca.proportion.of.variance)
    chart.title <- S4Vectors::metadata(TapestriExperiment)$pca.assay
  } else {
    knee.df <- data.frame("prop.variance" = S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.proportion.of.variance)
    chart.title <- S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$pca.assay
  }

  knee.df$cumulative.prop <- cumsum(knee.df$prop.variance)
  knee.df$index <- seq_len(nrow(knee.df))

  if (n.pcs > nrow(knee.df)) {
    warning("n.pcs exceeds number of principal components.")
    n.pcs <- nrow(knee.df)
  }

  knee.df <- knee.df[seq_len(n.pcs), ]
  knee.df$index <- as.factor(knee.df$index)

  g1 <- ggplot(knee.df) +
    geom_line(aes(x = .data$index, y = .data$prop.variance, group = 1)) +
    geom_point(aes(x = .data$index, y = .data$prop.variance)) +
    geom_col(aes(x = .data$index, y = .data$cumulative.prop), alpha = 0.3, width = 0.5) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 1)) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(
      x = "Principal Component",
      y = "Percent/cumulative Variance Explained",
      title = "Variance Explained by Principal Components",
      subtitle = chart.title
    )

  return(g1)
}

#' Cluster matrix data by UMAP
#'
#' Analyzes matrix data by UMAP and saves results to `reducedDims` slot of `TapestriObject`.
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param alt.exp Character, `altExp` to use, `NULL` uses top-level/main experiment. Default "alleleFrequency".
#' @param assay Character, `assay` to use. `NULL` (default) uses first-indexed assay. Not used when `use.pca.dims = TRUE`.
#' @param use.pca.dims Logical, if `TRUE`, uses experiment PCA, otherwise uses `assay` data. Default `TRUE`.
#' @param pca.dims Numeric, indices of PCs to use in UMAP. Default `NULL`.
#' @param ... Additional parameters to pass to [umap::umap()], e.g. for configuration (see [umap::umap.defaults()]).
#'
#' @return `TapestriExperiment` with UMAP embeddings saved to `reducedDims` slot of `altExp.`
#' @export
#'
#' @concept dimensional reduction
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' tap.object <- runPCA(tap.object, alt.exp = "alleleFrequency")
#' tap.object <- runUMAP(tap.object, pca.dims = 1:3)
runUMAP <- function(TapestriExperiment,
                    alt.exp = "alleleFrequency",
                    assay = NULL,
                    use.pca.dims = TRUE,
                    pca.dims = NULL, ...) {
  assay <- .SelectAssay(TapestriExperiment, alt.exp = alt.exp, assay = assay)

  if (alt.exp == assay) {
    cli::cli_alert_info("Running UMAP on: {.emph {alt.exp}}.")
  } else {
    cli::cli_alert_info("Running UMAP on: {.emph {alt.exp}} {.emph {assay}}.")
  }

  # PCA
  if (use.pca.dims) {
    if (is.null(alt.exp)) {
      umap.assay <- SingleCellExperiment::reducedDim(TapestriExperiment,
        type = "PCA"
      )
    } else {
      umap.assay <- SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(TapestriExperiment, alt.exp), type = "PCA")
    }

    if (is.null(pca.dims)) {
      cli::cli_abort("{.var pca.dims} not set.")
    } else {
      umap.assay <- umap.assay[, pca.dims]
    }
  } else {
    if (is.null(alt.exp)) {
      umap.assay <- SummarizedExperiment::assay(TapestriExperiment, assay)
    } else {
      umap.assay <- SummarizedExperiment::assay(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))
    }
  }

  # run UMAP
  umap.result <- umap::umap(umap.assay, ... = ...)
  umap.embeddings <- as.data.frame(umap.result$layout)
  colnames(umap.embeddings) <- c("umap.1", "umap.2")

  # save UMAP to main experiment
  if (is.null(alt.exp)) {
    SingleCellExperiment::reducedDim(TapestriExperiment, "UMAP", withDimnames = TRUE) <- umap.embeddings
    S4Vectors::metadata(TapestriExperiment)$"umap.assay" <-
      ifelse(use.pca.dims, "PCA", assay)
  } else { # save UMAP to alt experiment
    SingleCellExperiment::reducedDim(altExp(TapestriExperiment, alt.exp), "UMAP", withDimnames = TRUE) <- umap.embeddings
    S4Vectors::metadata(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))$"umap.assay" <- ifelse(use.pca.dims, "PCA", assay)
  }

  return(TapestriExperiment)
}

#' Scatter plot for dimensional reduction results
#'
#' Plots a scatter plot of the indicated dimensional reduction results.
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param dim.reduction Character, dimension reduction to plot, either "PCA" or "UMAP".
#' @param group.label Character, `colData` column for grouping samples by color. Default `NULL`.
#' @param alt.exp Character, `altExp` to use, `NULL` uses top-level/main experiment. Default "alleleFrequency".
#' @param dim.x Numeric, index of dimensional reduction data to plot on X axis. Default 1.
#' @param dim.y Numeric, index of dimensional reduction data to plot on Y axis. Default 2.
#'
#' @return ggplot2 object, scatter plot
#' @export
#'
#' @concept dimensional reduction
#' @concept plots
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' tap.object <- runPCA(tap.object, alt.exp = "alleleFrequency")
#' reducedDimPlot(tap.object, dim.reduction = "pca")
reducedDimPlot <- function(TapestriExperiment,
                           alt.exp = "alleleFrequency",
                           dim.reduction,
                           dim.x = 1,
                           dim.y = 2,
                           group.label = NULL) {
  .SelectAssay(TapestriExperiment, alt.exp = alt.exp)

  dim.reduction <- toupper(dim.reduction)

  if (is.null(alt.exp)) {
    to.plot <- SingleCellExperiment::reducedDim(TapestriExperiment, dim.reduction)
    to.plot <- to.plot[, c(dim.x, dim.y)]
  } else {
    to.plot <- SingleCellExperiment::reducedDim(altExp(TapestriExperiment, alt.exp), dim.reduction)
    to.plot <- to.plot[, c(dim.x, dim.y)]
  }

  if (is.null(group.label)) {
    group.label.data <- NULL
  } else {
    # check for group label
    if (!group.label %in% colnames(SummarizedExperiment::colData(TapestriExperiment))) {
      cli::cli_abort("group.label not found in colData.")
    }
    group.label.data <- SummarizedExperiment::colData(TapestriExperiment)[, group.label]
  }

  g1 <- simpleScatterPlot(
    x = to.plot[, 1],
    y = to.plot[, 2],
    labs.title = dim.reduction,
    labs.x = colnames(to.plot)[1],
    labs.y = colnames(to.plot[2]),
    group.label = group.label.data,
    group.label.legend = group.label
  )

  return(g1)
}

#' Cluster 2D data
#'
#' Clusters data using dbscan method and saves cluster assignments for each cell barcode to `colData`.
#' Generally used to assign clusters to UMAP projection after PCA and UMAP dimensional reduction.
#'
#' @param TapestriExperiment `TapestriExperiment` object
#' @param alt.exp Character, `altExp` slot to use. `NULL` uses top-level/main experiment. Default "alleleFrequency".
#' @param dim.reduction Character, reduced dimension data to use. Default "UMAP".
#' @param eps Numeric, `dbscan` `eps` parameter. Lower to increase cluster granularity. See [dbscan::dbscan()]. Default 0.8.
#' @param ... Additional parameters to pass to [dbscan::dbscan()].
#' @param dim.1 Numeric, index of data dimension to use. Default 1.
#' @param dim.2 Numeric, index of data dimension to use. Default 2.
#'
#' @return `TapestriExperiment` object with updated `colData` containing cluster assignments.
#' @export
#'
#' @concept dimensional reduction
#'
#' @seealso [dbscan::dbscan()]
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' tap.object <- runPCA(tap.object, alt.exp = "alleleFrequency")
#' tap.object <- runUMAP(tap.object, pca.dims = 1:3)
#' tap.object <- runClustering(tap.object, dim.reduction = "UMAP", eps = 0.8)
runClustering <- function(TapestriExperiment,
                          alt.exp = "alleleFrequency",
                          dim.reduction = "UMAP",
                          eps = 0.8,
                          dim.1 = 1,
                          dim.2 = 2,
                          ...) {
  dim.reduction <- toupper(dim.reduction)

  .SelectAssay(TapestriExperiment, alt.exp = alt.exp)

  cli::cli_alert_info("Finding clusters in: {.emph {alt.exp}} {.emph {dim.reduction}}")

  if (is.null(alt.exp)) {
    dbscan.assay <- SingleCellExperiment::reducedDim(TapestriExperiment, dim.reduction)
    dbscan.assay <- dbscan.assay[, c(dim.1, dim.2)]
  } else {
    dbscan.assay <- SingleCellExperiment::reducedDim(altExp(TapestriExperiment, alt.exp), dim.reduction)
    dbscan.assay <- dbscan.assay[, c(dim.1, dim.2)]
  }

  dbscan.result <- dbscan::dbscan(dbscan.assay, eps = eps, ...)
  dbscan.result.clusters <- data.frame(
    cell.barcode = rownames(dbscan.assay),
    cluster = as.factor(dbscan.result$cluster)
  )

  dbscan.result.clusters$cluster <- .FactorNumericalOrder(dbscan.result.clusters$cluster)

  ##
  # get and merge colData in main
  existing.cell.data <- as.data.frame(SummarizedExperiment::colData(TapestriExperiment))

  # drop existing clusters if they exist to allow overwriting
  existing.cell.data <- existing.cell.data[, which(colnames(existing.cell.data) != "cluster"), drop = FALSE]

  # merge result and existing colData
  updated.cell.data <- merge(existing.cell.data,
    dbscan.result.clusters,
    by = "cell.barcode",
    all.x = TRUE,
    sort = FALSE
  )

  # reorder to match colData
  rownames(updated.cell.data) <- updated.cell.data$cell.barcode
  updated.cell.data <- updated.cell.data[rownames(SingleCellExperiment::colData(TapestriExperiment)), ]

  # update TapestriExperiment
  SummarizedExperiment::colData(TapestriExperiment) <-
    S4Vectors::DataFrame(updated.cell.data)

  ##
  # get and merge colData in altExp
  existing.cell.data <- as.data.frame(SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, alt.exp)))

  # drop existing clusters if they exist to allow overwriting
  existing.cell.data <- existing.cell.data[, which(colnames(existing.cell.data) != "cluster"), drop = FALSE]

  # merge result and existing colData
  updated.cell.data <- merge(existing.cell.data,
    dbscan.result.clusters,
    by = "cell.barcode",
    all.x = TRUE,
    sort = FALSE
  )

  # reorder to match colData
  rownames(updated.cell.data) <- updated.cell.data$cell.barcode
  updated.cell.data <- updated.cell.data[rownames(SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, alt.exp))), ]

  # update TapestriExperiment
  SummarizedExperiment::colData(SingleCellExperiment::altExp(TapestriExperiment, alt.exp)) <- S4Vectors::DataFrame(updated.cell.data)

  return(TapestriExperiment)
}

.FactorNumericalOrder <- function(f) {
  f <- stats::reorder(f, f, FUN = length)
  f <- factor(f, levels = rev(levels(f)))
  levels(f) <- seq_len(length(levels(f)))

  return(f)
}
