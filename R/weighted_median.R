#' Smooth copy number values across chromosomes and chromosome arms
#'
#' `calcSmoothCopyNumber()` takes `copyNumber` slot values for probes on a chromosome and smooths them by median (default) for each chromosome
#' and chromosome arm, resulting in one copy number value per chromosome and chromosome arm for each cell barcode.
#' Cell-chromosome values are then discretized into integers by conventional rounding (1.5 <= x < 2.5 rounds to 2).
#' Smoothed copy number and discretized smoothed copy number values are stored as `smoothedCopyNumber` and `discreteCopyNumber` assays,
#' in `altExp` slots `smoothedCopyNumberByChr` for chromosome-level smoothing, and `smoothedCopyNumberByArm` for chromosome arm-level smoothing.
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
calcSmoothCopyNumber_FA <- function(TapestriExperiment, method = "median", control.copy.number = NULL, 
  sample.feature = "cluster", weight.range = c(0.5, 0.5)) {
  method <- tolower(method)

  if (method == "median") {
    smooth.func <- stats::median
    S4Vectors::metadata(TapestriExperiment)$smoothing.method <- "median"
  } else if (method == "mean") {
    smooth.func <- mean
    S4Vectors::metadata(TapestriExperiment)$smoothing.method <- "mean"
  } else if (method == "weighted.median") {
      smooth.func <- NULL
      S4Vectors::metadata(TapestriExperiment)$smoothing.method <- "weighted.median"
      S4Vectors::metadata(TapestriExperiment)$smoothing.weights <- weight.range
  } else {
    cli::cli_abort("{.var method} {.q {method}}, not recognized. Please use {.q mean} or {.q median}.")
  }
  
  cli::cli_progress_step("Smoothing copy number by {method}...", )
  
  ploidy.counts <- SummarizedExperiment::assay(TapestriExperiment, "copyNumber")
  
  tap.exp.row.data <- as.data.frame(SummarizedExperiment::rowData(TapestriExperiment))[, c("probe.id", "chr", "arm", "cytoband")]
  tap.exp.row.data$cytoband <- paste0('chr',tap.exp.row.data$chr,tap.exp.row.data$cytoband)


  #get weight for each probe
  if(method == "weighted.median"){

      tap.exp.row.data <- tap.exp.row.data %>% dplyr::left_join(control.copy.number, by = "cytoband")
        
      tap.exp.row.data$probe.weight <- tap.exp.row.data[, c("probe.id", "copy.number", "sample.label")] %>% 
          purrr::pmap(function(probe.id, copy.number, sample.label){
              barcodes <- SingleCellExperiment::colData(TapestriExperiment)[SingleCellExperiment::colData(TapestriExperiment)[,{{sample.feature}}] == sample.label,"cell.barcode"]
              copy.number.values <- SummarizedExperiment::assay(TapestriExperiment, "copyNumber")[probe.id,barcodes]
              lower.bound <- copy.number - weight.range[1]
              upper.bound <- copy.number + weight.range[2]
              accuracy <- round(sum(copy.number.values > lower.bound & copy.number.values < upper.bound)/length(copy.number.values), 3)
              return(accuracy)
          }) %>% unlist()
      
      S4Vectors::metadata(TapestriExperiment)$probe.weights <- tap.exp.row.data
  }

  ploidy.tidy <- ploidy.counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("probe.id") %>%
    tidyr::pivot_longer(
      cols = !tidyr::matches("probe.id"),
      names_to = "cell.barcode",
      values_to = "ploidy"
    ) 


  if(method == "median" | method == "mean"){

    ploidy.tidy <- ploidy.tidy %>% dplyr::left_join(tap.exp.row.data[, c("probe.id", "chr", "arm", "cytoband")], by = "probe.id")

    # whole-chromosome 
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


    # arms
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


    # cytobands
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

  }else if(method == "weighted.median"){  #weighed median smoothing 
      
      ploidy.tidy <- ploidy.tidy %>% dplyr::left_join(tap.exp.row.data[, c("probe.id", "chr", "arm", "probe.weight")], by = "probe.id")
      
      # whole chromosome
      smoothed.ploidy.chr <- ploidy.tidy %>% 
          dplyr::group_split(.data$chr, .data$cell.barcode, .keep = TRUE) %>% 
          lapply(function(x){
              result <- round(matrixStats::weightedMedian(x = x$ploidy, w = x$probe.weight), 3)
              return(list(x$cell.barcode[1], x$chr[1], result))
          }) %>% 
          purrr::list_transpose()
      
      smoothed.ploidy.chr <- data.frame(cell.barcode = smoothed.ploidy.chr[[1]], 
                                        feature.id = smoothed.ploidy.chr[[2]], 
                                        value = as.numeric(smoothed.ploidy.chr[[3]])) 
      
      smoothed.ploidy.chr <- tidyr::pivot_wider(smoothed.ploidy.chr, names_from = .data$cell.barcode, values_from = .data$value) %>% 
          tibble::column_to_rownames("feature.id")
      
      # reorder to match input matrix
      smoothed.ploidy.chr <- smoothed.ploidy.chr[, colnames(ploidy.counts)]


      # arms
      smoothed.ploidy.arm <- ploidy.tidy %>% 
          dplyr::group_split(.data$arm, .data$cell.barcode, .keep = TRUE) %>% 
          lapply(function(x){
              result <- round(matrixStats::weightedMedian(x = x$ploidy, w = x$probe.weight), 3)
              return(list(x$cell.barcode[1], x$arm[1], result))
          }) %>% 
          purrr::list_transpose()
      
      smoothed.ploidy.arm <- data.frame(cell.barcode = smoothed.ploidy.arm[[1]], 
                                        feature.id = smoothed.ploidy.arm[[2]], 
                                        value = as.numeric(smoothed.ploidy.arm[[3]])) 
      
      smoothed.ploidy.arm <- tidyr::pivot_wider(smoothed.ploidy.arm, names_from = .data$cell.barcode, values_from = .data$value) %>% 
          tibble::column_to_rownames("feature.id")

      # reorder to match input matrix
      smoothed.ploidy.arm <- smoothed.ploidy.arm[, colnames(ploidy.counts)]


      # cytobands
      smoothed.ploidy.cytob <- ploidy.tidy %>% 
          dplyr::group_split(.data$cytobands, .data$cell.barcode, .keep = TRUE) %>% 
          lapply(function(x){
              result <- round(matrixStats::weightedMedian(x = x$ploidy, w = x$probe.weight), 3)
              return(list(x$cell.barcode[1], x$cytobands[1], result))
          }) %>% 
          purrr::list_transpose()
      
      smoothed.ploidy.cytob <- data.frame(cell.barcode = smoothed.ploidy.cytob[[1]], 
                                          feature.id = smoothed.ploidy.cytob[[2]], 
                                          value = as.numeric(smoothed.ploidy.cytob[[3]])) 
      
      smoothed.ploidy.cytob <- tidyr::pivot_wider(smoothed.ploidy.cytob, names_from = .data$cell.barcode, values_from = .data$value) %>% 
          tibble::column_to_rownames("feature.id")

      # reorder to march input matrix
      smoothed.ploidy.cytob <- smoothed.ploidy.cytob[, colnames(ploidy.counts)]
     
  }

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