#' Normalize raw counts
#'
#' Normalizes raw counts from `counts` slot in `TapestriExperiment` and returns the object with normalized counts in the `normcounts` slot.
#' Also calculates the standard deviation for each probe using normalized counts and adds it to `rowData`.
#'
#' "mb" method performs the same normalization scheme as in Mission Bio's mosaic package for python:
#' Counts for each barcode are normalized relative to their barcode's mean and probe counts are normalized relative to their probe's median.
#' "libNorm" method preforms library size normalization, returning the proportion of counts of each probe within a cell.
#' The proportion is multiplied by `scaling.factor` if provided.
#'
#' @param TapestriExperiment `TapestriExperiment` object.
#' @param method Character, normalization method. Default "mb".
#' @param scaling.factor Numeric, optional number to scale normalized counts if `method == "libNorm"`. Default `NULL`.
#'
#' @return `TapestriExperiment` object with normalized counts added to `normcounts` slot.
#' @export
#'
#' @concept copy number
#'
#' @importFrom stats median sd
#'
#' @examples
#' tap.object <- newTapestriExperimentExample() # example TapestriExperiment object
#' tap.object <- calcNormCounts(tap.object)
calcNormCounts <- function(TapestriExperiment,
                           method = "kt",
                           scaling.factor = 1000,
                           filter.bc = TRUE,
                           limits.filter.bc = c(0.5,2.5)) {
  method <- tolower(method)

  if(filter.bc){
    upper_limit <- median(TapestriExperiment$total.reads)+limits.filter.bc[2]*sd(TapestriExperiment$total.reads)
    lower_limit <- median(TapestriExperiment$total.reads)-limits.filter.bc[1]*sd(TapestriExperiment$total.reads)
    cell.bc <- colData(TapestriExperiment) %>% as_tibble() %>% 
      filter(total.reads <= upper_limit & total.reads >= lower_limit) %>% select(cell.barcode)
    
    TapestriExperiment <- TapestriExperiment[, colData(TapestriExperiment)$cell.barcode %in% cell.bc$cell.barcode]
  }
  
  raw.count.matrix <- SummarizedExperiment::assay(TapestriExperiment, "counts")

  if (any(is.na(raw.count.matrix))) {
    warning("NAs found in count data. normcounts may also contain NAs.")
  }

  if (method == "kt") {
    read.counts.normal <- .ktNormCounts(raw.count.matrix, scaling.factor)
  } else if (method == "kt2") {
    read.counts.normal <- .kt2NormCounts(raw.count.matrix)
  } else if (method == "mb") {
    read.counts.normal <- .MBNormCounts(raw.count.matrix)
  } else if (method == "libnorm") {
    read.counts.normal <- .LibSizeNorm(raw.count.matrix,
      scaling.factor = scaling.factor
    )
  } else {
    warning("Method not recognized. Set method to 'mb' or libNorm'")
  }

  # add to normalized counts slot
  normcounts(TapestriExperiment) <- read.counts.normal

  # calculate norm count SD for each amplicon and add to rowData
  current.probe.data <- SummarizedExperiment::rowData(TapestriExperiment)
  current.probe.order <- rownames(current.probe.data)
  new.probe.data <- apply(read.counts.normal, 1, sd)
  names(new.probe.data) <- rownames(read.counts.normal)
  new.probe.data <- new.probe.data[current.probe.order]
  SummarizedExperiment::rowData(TapestriExperiment)$norm.count.sd <- new.probe.data

  return(TapestriExperiment)
}

.ktNormCounts <- function(input.matrix, scaling.factor){
  # lib size normalization 
  input.matrix <- apply(input.matrix, 2, function(x)(x)/sum(x)) * scaling.factor
  # probe normalization
  matrix.normal <- apply(input.matrix, 1, function(x)(x+1)/sum(x)) * scaling.factor
  matrix.normal <- t(matrix.normal)

  return(matrix.normal)
}

.kt2NormCounts <- function(input.matrix){
  # lib size normalization 
  input.matrix <- apply(input.matrix, 2, function(x)(x)/sum(x)) 
  input.matrix <- input.matrix * (1/median(input.matrix[input.matrix != 0])) # estimate scaling.factor
  # probe normalization
  matrix.normal <- apply(input.matrix, 1, function(x)(x+1)/median(x))
  matrix.normal <- t(matrix.normal)
  
  return(matrix.normal)
}

.MBNormCounts <- function(input.matrix) {
  # get "good barcodes", barcodes that have at least 10% the counts of the 11th barcode ranked for highest number of counts
  barcode.sums <- apply(input.matrix, MARGIN = 2, sum)
  barcode.sorted <- sort(barcode.sums, decreasing = TRUE)
  good.barcodes <- barcode.sums > (barcode.sorted[11] / 10)

  # normalize barcodes relative to barcode means
  barcode.means <- apply(input.matrix, 2, mean) + 1
  matrix.normal <- sweep(input.matrix, 2, barcode.means, "/")

  # normalize probes relative to probe median. medians calculated using "good barcodes"
  probe.medians <- apply(matrix.normal[, good.barcodes], 1, median) + 0.05
  matrix.normal <- sweep(matrix.normal, 1, probe.medians, "/")

  # scale all counts by 2X for diploid baseline
  matrix.normal <- matrix.normal * 2

  return(matrix.normal)
}

.LibSizeNorm <- function(input.matrix, scaling.factor = NULL) {
  column.sums <- colSums(input.matrix)
  matrix.normal <- sweep(input.matrix, 2, column.sums, "/")

  if (!is.null(scaling.factor)) {
    matrix.normal <- matrix.normal * scaling.factor
  }

  return(matrix.normal)
}
