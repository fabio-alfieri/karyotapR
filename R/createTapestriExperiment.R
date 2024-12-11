#' Create `TapestriExperiment` object from Tapestri Pipeline output
#'
#' `createTapestriExperiment()` constructs a `TapestriExperiment` container object from data stored in the `.h5` file output by the Tapestri Pipeline.
#' Read count matrix (probe x cell barcode) is stored in the "counts" `assay` slot of the top-level experiment.
#' Allele frequency matrix (variant x cell barcode) is stored in the "alleleFrequency" `assay` slot of the "alleleFrequency" `altExp` (alternative experiment) slot.
#' `panel.id` is an optional shortcut to set special probe identities for specific custom panels.
#'
#' # Panel ID Shortcuts
#' `panel.id` is an optional shortcut to set the `barcodeProbe` and `grnaProbe` slots in `TapestriExperiment` for specific custom Tapestri panels.
#' ## CO261
#' - `barcodeProbe` = "not specified"
#' - `grnaProbe` = "not specified"
#'
#' ## CO293
#' - `barcodeProbe` = "AMPL205334"
#' - `grnaProbe` = "AMPL205666"
#'
#' ## CO610
#' - `barcodeProbe` = "CO610_AMP351"
#' - `grnaProbe` = "CO610_AMP350"
#' 
#' ## CO810
#' - `barcodeProbe` = "TAMPL46684"
#' - `grnaProbe` = "TAMPL46683"
#'
#' # Automatic Operations
#' ## Raw Data
#' Read count and allele frequency matrices are imported to their appropriate slots as described above.
#' `filter.variants == TRUE` (default) only loads allele frequency variants that have passed internal filters in the Tapestri Pipeline.
#' This greatly reduces the number of variants from tens of thousands to hundreds of likely more consequential variants,
#' saving RAM and reducing operation time.
#'
#' ## Metadata
#' Several metadata sets are copied or generated and then stored in the appropriate `TapestriExperiment` slot during construction.
#' - Probe panel metadata stored in the `.h5` file are copied to `rowData`.
#' - Basic QC stats (e.g. total number of reads per probe) are added to `rowData.`
#' - Basic QC stats (e.g. total number of reads per cell barcode) are added to `colData.`
#' - Experiment-level metadata is stored in `metadata`.
#'
#' # Optional Operations
#' Two additional major operations are called by default during `TapestriExperiment` construction for convenience.
#' `get.cytobands == TRUE` (default) calls [getCytobands()], which retrieves the chromosome arm and cytoband for each probe based on stored positional data and saves them in `rowData`.
#' Some downstream smoothing and plotting functions may fail if chromosome arms are not present in `rowData`, so this generally should always be run.
#' `move.non.genome.probes` calls [moveNonGenomeProbes()], which moves probes corresponding to the specified tags to `altExp` (alternative experiment) slots in the `TapestriExperiment` object.
#' The exception is probes on chromosome Y; CNVs of chrY are more rare, so we move it to an `altExp` for separate analysis.
#' Probes corresponding to the `barcodeProbe` and `grnaProbe` slots, which are specified by the `panel.id` shortcut or manually (see [Custom Slot Getters and Setters]),
#' are automatically moved to `altExp` by this operation as well.
#' If such probes are not present, the function will only generate a warning message, so it is always safe (and recommended) to run by default.
#' Any remaining probes that are not targeting a human chromosome and are not specified by the shortcut tags are moved to the `otherProbeCounts` slot.
#'
#' @param h5.filename File path for `.h5` file from Tapestri Pipeline output.
#' @param panel.id Character, Tapestri panel ID, either CO261, CO293, CO610, CO810, or `NULL`. Initializes `barcodeProbe` and `grnaProbe` slots. Default `NULL`.
#' @param get.cytobands Logical, if `TRUE` (default), retrieve and add chromosome cytobands and chromosome arms to `rowData` (probe metadata).
#' @param genome Character, reference genome for pulling cytoband coordinates and chromosome arm labels (see [getCytobands()]). Only "hg19" (default) is currently supported.
#' @param move.non.genome.probes Logical, if `TRUE` (default), move counts and metadata from non-genomic probes to `altExp` slots (see [moveNonGenomeProbes()]).
#' @param filter.variants Logical, if `TRUE` (default), only stores variants that have passed Tapestri Pipeline filters.
#' @param verbose Logical, if `TRUE` (default), metadata is output in message text.
#'
#' @return `TapestriExperiment` object containing data from Tapestri Pipeline output.
#' @export
#'
#' @import SingleCellExperiment
#'
#' @seealso [moveNonGenomeProbes()], [getCytobands()], which are run as part of this function by default.
#'
#' @concept build experiment
#'
#' @examples
#' \dontrun{
#' tapExperiment <- createTapestriExperiment("myh5file.h5", "CO293")
#' }
createTapestriExperiment <- function(h5.filename,
                                     panel.id = NULL,
                                     get.cytobands = TRUE,
                                     genome = "hg19",
                                     move.non.genome.probes = TRUE,
                                     filter.variants = TRUE,
                                     verbose = TRUE) {

    # read panel ID
  panel.id.output <- .GetPanelID(panel.id = panel.id)
  barcodeProbe <- panel.id.output[["barcode.probe"]]
  grnaProbe <- panel.id.output[["grna.probe"]]


  # import data
  tapestri.h5 <- rhdf5::H5Fopen(file.path(h5.filename))

  # construct object
  read.counts.raw <- t(matrix(
    data = tapestri.h5$"/assays/dna_read_counts/layers/read_counts",
    ncol = length(tapestri.h5$"/assays/dna_read_counts/ca/id"),
    byrow = TRUE
  ))

  read.counts.raw.colData <- S4Vectors::DataFrame(
    cell.barcode = as.character(tapestri.h5$"/assays/dna_read_counts/ra/barcode"),
    row.names = as.character(tapestri.h5$"/assays/dna_read_counts/ra/barcode")
  )

  read.counts.raw.rowData <- S4Vectors::DataFrame(
    probe.id = as.character(tapestri.h5$"/assays/dna_read_counts/ca/id"),
    chr = tapestri.h5$"/assays/dna_read_counts/ca/CHROM",
    start.pos = as.numeric(tapestri.h5$"/assays/dna_read_counts/ca/start_pos"),
    end.pos = as.numeric(tapestri.h5$"/assays/dna_read_counts/ca/end_pos"),
    row.names = as.character(tapestri.h5$"/assays/dna_read_counts/ca/id")
  )

  chr.order <- getChrOrder(read.counts.raw.rowData$chr) # reorder amplicon metadata by chromosome order

  read.counts.raw.rowData <- read.counts.raw.rowData[chr.order, ]

  read.counts.raw <- read.counts.raw[chr.order, ]

  read.counts.raw.rowData$chr <- factor(read.counts.raw.rowData$chr,
                                        levels = unique(read.counts.raw.rowData$chr))

  # basic metrics
  read.counts.raw.colData$total.reads <- colSums(read.counts.raw) # total reads per cell
  read.counts.raw.rowData$total.reads <- rowSums(read.counts.raw) # total reads per probe
  read.counts.raw.rowData$median.reads <- apply(read.counts.raw, 1, median) # median reads per probe
  mean.reads <- round(mean(colMeans(read.counts.raw)), 2) # mean reads/cell/probe(or amplicon)

  # report data
  if (verbose) {
    cli::cli_h1("Loading Tapestri Experiment")
    cli::cli_bullets(c(
      "*" = 'Sample Name: {tapestri.h5$"/assays/dna_read_counts/metadata/sample_name"}',
      "*" = 'Pipeline Panel Name: {tapestri.h5$"/assays/dna_read_counts/metadata/panel_name"}',
      "*" = 'Pipeline Version: {tapestri.h5$"/assays/dna_read_counts/metadata/pipeline_version"}',
      "*" = 'Date Created: {tapestri.h5$"/metadata/date_created"}'
    ))
    cli::cli_h2("Metrics")
    cli::cli_bullets(c(
      "*" = 'Number of Cells: {tapestri.h5$"/assays/dna_read_counts/metadata/n_cells"}',
      "*" = 'Number of Probes: {tapestri.h5$"/assays/dna_read_counts/metadata/n_amplicons"}',
      "*" = "Mean Reads per Cell per Probe: {mean.reads}"
    ))
  }

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
    colData = read.counts.raw.colData,
    rowData = read.counts.raw.rowData
  )

  tapestri.object <- .TapestriExperiment(sce)

  SingleCellExperiment::mainExpName(tapestri.object) <- "CNV"

  # save experiment metadata
  S4Vectors::metadata(tapestri.object)$sample.name <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/sample_name")
  S4Vectors::metadata(tapestri.object)$pipeline.panel.name <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/panel_name")
  S4Vectors::metadata(tapestri.object)$pipeline.version <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/pipeline_version")
  S4Vectors::metadata(tapestri.object)$number.of.cells <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/n_cells")
  S4Vectors::metadata(tapestri.object)$number.of.probes <- as.character(tapestri.h5$"/assays/dna_read_counts/metadata/n_amplicons")
  S4Vectors::metadata(tapestri.object)$date.h5.created <- as.character(tapestri.h5$"/metadata/date_created")
  S4Vectors::metadata(tapestri.object)$mean.reads.per.cell.per.probe <- as.character(mean.reads)

  # apply panel ID probe shortcuts
  tapestri.object@barcodeProbe <- barcodeProbe
  tapestri.object@grnaProbe <- grnaProbe

  # if keeping only variants that passed thresholds
  if (filter.variants == TRUE) {
    # variant filtering; filtered == TRUE in h5 object means "filtered out".
    filtered.variants <- as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered")
    filtered.variants <- !filtered.variants

    # variant allele frequency data
    variant.metadata <- data.frame(
      filtered = as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered")[filtered.variants],
      chr = as.character(tapestri.h5$"/assays/dna_variants/ca/CHROM")[filtered.variants],
      position = as.numeric(tapestri.h5$"/assays/dna_variants/ca/POS")[filtered.variants],
      quality = as.numeric(tapestri.h5$"/assays/dna_variants/ca/QUAL")[filtered.variants],
      amplicon.id = as.character(tapestri.h5$"/assays/dna_variants/ca/amplicon")[filtered.variants],
      variant.id = as.character(tapestri.h5$"/assays/dna_variants/ca/id")[filtered.variants],
      ado.rate = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_rate")[filtered.variants],
      reference.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/REF")[filtered.variants],
      alternate.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/ALT")[filtered.variants],
      ado.gt.cells = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_gt_cells")[filtered.variants]
    )

    af.matrix <- tapestri.h5$"/assays/dna_variants/layers/AF"[filtered.variants, ]
    dimnames(af.matrix) <- list(variant.metadata$variant.id, tapestri.h5$"/assays/dna_read_counts/ra/barcode")
  } else {
    # variant allele frequency data, no filtering, keep all variants
    variant.metadata <- data.frame(
      filtered = as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered"),
      chr = as.character(tapestri.h5$"/assays/dna_variants/ca/CHROM"),
      position = as.numeric(tapestri.h5$"/assays/dna_variants/ca/POS"),
      quality = as.numeric(tapestri.h5$"/assays/dna_variants/ca/QUAL"),
      amplicon.id = as.character(tapestri.h5$"/assays/dna_variants/ca/amplicon"),
      variant.id = as.character(tapestri.h5$"/assays/dna_variants/ca/id"),
      ado.rate = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_rate"),
      reference.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/REF"),
      alternate.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/ALT"),
      ado.gt.cells = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_gt_cells")
    )

    af.matrix <- tapestri.h5$"/assays/dna_variants/layers/AF"
    dimnames(af.matrix) <- list(variant.metadata$variant.id, tapestri.h5$"/assays/dna_read_counts/ra/barcode")
  }

  variant.metadata$chr <- factor(variant.metadata$chr, unique(variant.metadata$chr))
  variant.metadata$allelefreq.sd <- apply(af.matrix, 1, stats::sd)
  rownames(variant.metadata) <- variant.metadata$variant.id

  allele.frequency <- SingleCellExperiment::SingleCellExperiment(list(alleleFrequency = af.matrix),
    rowData = S4Vectors::DataFrame(variant.metadata),
    colData = S4Vectors::DataFrame(cell.barcode = colnames(af.matrix))
  )

  allele.frequency <- .TapestriExperiment(allele.frequency)
  allele.frequency@barcodeProbe <- barcodeProbe
  allele.frequency@grnaProbe <- grnaProbe

  SingleCellExperiment::altExp(tapestri.object, "alleleFrequency", withDimnames = TRUE) <- allele.frequency

  # close h5
  rhdf5::H5Fclose(tapestri.h5)

  # get cytobands
  if (get.cytobands) {
    tapestri.object <- getCytobands(tapestri.object, verbose = verbose)
  }

  # move non-genomic probes to altExp slots
  if (!identical(move.non.genome.probes, FALSE)) {
    tapestri.object <- moveNonGenomeProbes(tapestri.object)
  }

  return(tapestri.object)
}

.GetPanelID <- function(panel.id) {
  # read panel ID
  if (is.null(panel.id)) {
    barcodeProbe <- "not specified"
    grnaProbe <- "not specified"
  } else if (panel.id == "CO293") {
    barcodeProbe <- "AMPL205334"
    grnaProbe <- "AMPL205666"
  } else if (panel.id == "CO261") {
    barcodeProbe <- "not specified"
    grnaProbe <- "not specified"
  } else if (panel.id == "CO610") {
    barcodeProbe <- "CO610_AMP351"
    grnaProbe <- "CO610_AMP350"
  } else if (panel.id == "CO810") {
    barcodeProbe <- "TAMPL46684"
    grnaProbe <- "TAMPL46683"
  } else {
    cli::cli_abort("{.var panel.id} {.q {panel.id}} is not recognized. Please specify CO261, CO293, CO610, CO810, or NULL for manual settings.", )
  }

  return(list(barcode.probe = barcodeProbe, grna.probe = grnaProbe))
}

# create TapestriExperiment with manual counts and existing metadata for simulations
.createTapestriExperiment.sim <- function(counts = NULL,
                                          probe.metadata = NULL,
                                          genome = "hg19") {
  if (is.null(counts)) {
    cli::cli_abort("No counts matrix supplied.")
  }

  if (is.null(dimnames(counts))) {
    cli::cli_abort("No dimnames on counts matrix supplied.")
  }

  # read counts
  read.counts.raw <- as.matrix(counts)

  read.counts.raw.colData <- S4Vectors::DataFrame(
    cell.barcode = colnames(read.counts.raw),
    row.names = colnames(read.counts.raw)
  )

  # panel.metadata
  probe.metadata <- probe.metadata[, c("probe.id", "chr", "start.pos", "end.pos")]

  if (!all(rownames(read.counts.raw) == probe.metadata$probe.id)) {
    cli::cli_abort("Probe IDs in metadata do not match rownames of counts.")
  }

  read.counts.raw.rowData <- S4Vectors::DataFrame(
    probe.id = probe.metadata$probe.id,
    chr = probe.metadata$chr,
    start.pos = probe.metadata$start.pos,
    end.pos = probe.metadata$end.pos,
    row.names = rownames(probe.metadata)
  )

  # basic metrics
  read.counts.raw.colData$total.reads <- colSums(read.counts.raw) # total reads per cell
  read.counts.raw.rowData$total.reads <- rowSums(read.counts.raw) # total reads per probe
  read.counts.raw.rowData$median.reads <- apply(read.counts.raw, 1, median) # median reads per probe
  mean.reads <- round(mean(colMeans(read.counts.raw)), 2) # mean reads/cell/probe(or amplicon)

  ## BUILD OBJECT

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
    colData = read.counts.raw.colData,
    rowData = read.counts.raw.rowData
  )

  tapestri.object <- .TapestriExperiment(sce)

  SingleCellExperiment::mainExpName(tapestri.object) <- "CNV"

  # save experiment metadata
  S4Vectors::metadata(tapestri.object)$sample.name <- ""
  S4Vectors::metadata(tapestri.object)$pipeline.panel.name <- ""
  S4Vectors::metadata(tapestri.object)$pipeline.version <- ""
  S4Vectors::metadata(tapestri.object)$number.of.cells <- ncol(read.counts.raw)
  S4Vectors::metadata(tapestri.object)$number.of.probes <- nrow(read.counts.raw)
  S4Vectors::metadata(tapestri.object)$date.h5.created <- ""
  S4Vectors::metadata(tapestri.object)$mean.reads.per.cell.per.probe <- as.character(mean.reads)

  tapestri.object <- suppressMessages(getCytobands(tapestri.object, genome = genome))

  return(tapestri.object)
}

# create TapestriExperiment with manual counts
.createTapestriExperiment.manual <- function(counts = NULL, af.matrix = NULL, panel.id = NULL, get.cytobands = TRUE, genome = "hg19", move.non.genome.probes = TRUE,
                                             verbose = TRUE) {
  if (is.null(counts)) {
    cli::cli_abort("No counts matrix supplied.")
  }

  if (is.null(dimnames(counts))) {
    cli::cli_abort("No dimnames on counts matrix supplied.")
  }

  # read panel ID
  panel.id.output <- .GetPanelID(panel.id = panel.id)
  barcodeProbe <- panel.id.output[["barcode.probe"]]
  grnaProbe <- panel.id.output[["grna.probe"]]

  # read counts
  read.counts.raw <- as.matrix(counts)

  read.counts.raw.colData <- S4Vectors::DataFrame(
    cell.barcode = colnames(read.counts.raw),
    row.names = colnames(read.counts.raw)
  )

  # subset probe metadata in case not all probes present in count data

  if (panel.id == "CO261") {
    probe.metadata <- co261.metadata
  } else if (panel.id == "CO293") {
    probe.metadata <- co293.metadata
  } else if (panel.id == "CO610") {
    probe.metadata <- co610.metadata
  } else if (panel.id == "CO810") {
      probe.metadata <- co810.metadata
  }

  if (!all(rownames(read.counts.raw) %in% probe.metadata$probe.id)) {
    cli::cli_abort("Unknown probe ID in read counts.")
  }

  probe.metadata <- probe.metadata[probe.metadata$probe.id %in% rownames(read.counts.raw), ]

  read.counts.raw.rowData <- S4Vectors::DataFrame(
    probe.id = probe.metadata$probe.id,
    chr = probe.metadata$chr,
    start.pos = probe.metadata$start.pos,
    end.pos = probe.metadata$end.pos,
    row.names = rownames(probe.metadata)
  )

  chr.order <- getChrOrder(read.counts.raw.rowData$chr) # reorder amplicon metadata by chromosome order

  read.counts.raw.rowData <- read.counts.raw.rowData[chr.order, ]

  read.counts.raw <- read.counts.raw[chr.order, ]

  read.counts.raw.rowData$chr <- factor(read.counts.raw.rowData$chr, levels = unique(read.counts.raw.rowData$chr))

  # basic metrics
  read.counts.raw.colData$total.reads <- colSums(read.counts.raw) # total reads per cell
  read.counts.raw.rowData$total.reads <- rowSums(read.counts.raw) # total reads per probe
  read.counts.raw.rowData$median.reads <- apply(read.counts.raw, 1, median) # median reads per probe
  mean.reads <- round(mean(colMeans(read.counts.raw)), 2) # mean reads/cell/probe(or amplicon)

  ## BUILD OBJECT

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
    colData = read.counts.raw.colData,
    rowData = read.counts.raw.rowData
  )

  tapestri.object <- .TapestriExperiment(sce)

  SingleCellExperiment::mainExpName(tapestri.object) <- "CNV"

  # save experiment metadata
  S4Vectors::metadata(tapestri.object)$sample.name <- ""
  S4Vectors::metadata(tapestri.object)$pipeline.panel.name <- panel.id
  S4Vectors::metadata(tapestri.object)$pipeline.version <- ""
  S4Vectors::metadata(tapestri.object)$number.of.cells <- ncol(read.counts.raw)
  S4Vectors::metadata(tapestri.object)$number.of.probes <- nrow(read.counts.raw)
  S4Vectors::metadata(tapestri.object)$date.h5.created <- ""
  S4Vectors::metadata(tapestri.object)$mean.reads.per.cell.per.probe <- as.character(mean.reads)

  # apply panel ID probe shortcuts
  tapestri.object@barcodeProbe <- barcodeProbe
  tapestri.object@grnaProbe <- grnaProbe

  # get cytobands
  if (get.cytobands) {
    tapestri.object <- getCytobands(tapestri.object)
  }

  # move non-genomic probes to altExp slots
  if (!identical(move.non.genome.probes, FALSE)) {
    tapestri.object <- moveNonGenomeProbes(tapestri.object)
  }

  # report data
  if (verbose) {
    cli::cli_h1("Loading Tapestri Experiment")
    cli::cli_bullets(c(
      "*" = "Sample Name: {S4Vectors::metadata(tapestri.object)$sample.name}",
      "*" = "Pipeline Panel Name: {S4Vectors::metadata(tapestri.object)$pipeline.panel.name}",
      "*" = "Pipeline Version: {S4Vectors::metadata(tapestri.object)$pipeline.version}",
      "*" = "Date Created: {S4Vectors::metadata(tapestri.object)$date.h5.created}"
    ))
    cli::cli_h2("Metrics")
    cli::cli_bullets(c(
      "*" = "Number of Cells: {S4Vectors::metadata(tapestri.object)$number.of.cells}",
      "*" = "Number of Probes: {S4Vectors::metadata(tapestri.object)$number.of.probes}",
      "*" = "Mean Reads per Cell per Probe: {S4Vectors::metadata(tapestri.object)$mean.reads.per.cell.per.probe}"
    ))
  }

  return(tapestri.object)

  # allele frequency

  # not supported yet

  #     # variant filtering; filtered == TRUE in h5 object means "filtered out".
  #     filtered.variants <- as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered")
  #     filtered.variants <- !filtered.variants
  #
  #     # variant allele frequency data
  #     variant.metadata <- data.frame(
  #         filtered = as.logical(tapestri.h5$"/assays/dna_variants/ca/filtered")[filtered.variants],
  #         chr = as.character(tapestri.h5$"/assays/dna_variants/ca/CHROM")[filtered.variants],
  #         position = as.numeric(tapestri.h5$"/assays/dna_variants/ca/POS")[filtered.variants],
  #         quality = as.numeric(tapestri.h5$"/assays/dna_variants/ca/QUAL")[filtered.variants],
  #         amplicon.id = as.character(tapestri.h5$"/assays/dna_variants/ca/amplicon")[filtered.variants],
  #         variant.id = as.character(tapestri.h5$"/assays/dna_variants/ca/id")[filtered.variants],
  #         ado.rate = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_rate")[filtered.variants],
  #         reference.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/REF")[filtered.variants],
  #         alternate.allele = as.character(tapestri.h5$"/assays/dna_variants/ca/ALT")[filtered.variants],
  #         ado.gt.cells = as.numeric(tapestri.h5$"/assays/dna_variants/ca/ado_gt_cells")[filtered.variants]
  #     )
  #
  #     af.matrix <- tapestri.h5$"/assays/dna_variants/layers/AF"[filtered.variants, ]
  #     dimnames(af.matrix) <- list(variant.metadata$variant.id, tapestri.h5$"/assays/dna_read_counts/ra/barcode")
  #
  # variant.metadata$chr <- factor(variant.metadata$chr, unique(variant.metadata$chr))
  # variant.metadata$allelefreq.sd <- apply(af.matrix, 1, stats::sd)
  # rownames(variant.metadata) <- variant.metadata$variant.id
  #
  # allele.frequency <- SingleCellExperiment::SingleCellExperiment(list(alleleFrequency = af.matrix),
  #                                                                rowData = S4Vectors::DataFrame(variant.metadata),
  #                                                                colData = S4Vectors::DataFrame(cell.barcode = colnames(af.matrix))
  # )
  #
  # allele.frequency <- .TapestriExperiment(allele.frequency)
  # allele.frequency@barcodeProbe <- barcodeProbe
  # allele.frequency@grnaProbe <- grnaProbe
  #
  # SingleCellExperiment::altExp(tapestri.object, "alleleFrequency", withDimnames = TRUE) <- allele.frequency

  return(tapestri.object)
}

# create Dummy TapestriExperiment with CO293 metadata
newDummyTapestriExperiment <- function() {
  panel.id <- "CO293"

  # read panel ID
  panel.id.output <- .GetPanelID(panel.id = panel.id)
  barcodeProbe <- panel.id.output[["barcode.probe"]]
  grnaProbe <- panel.id.output[["grna.probe"]]

  # raw counts, 300 cells, panel CO293
  read.counts.raw <- matrix(data = 1:91200, ncol = 300)

  read.counts.raw.colData <- S4Vectors::DataFrame(
    cell.barcode = as.character(paste0("cell_", 1:300)),
    row.names = as.character(paste0("cell_", 1:300))
  )

  read.counts.raw.rowData <- S4Vectors::DataFrame(
    probe.id = as.character(co293.metadata$probe.id),
    chr = co293.metadata$chr,
    start.pos = as.numeric(co293.metadata$start.pos),
    end.pos = as.numeric(co293.metadata$end.pos),
    row.names = as.character(co293.metadata$probe.id)
  )

  chr.order <- getChrOrder(read.counts.raw.rowData$chr) # reorder amplicon metadata by chromosome order

  read.counts.raw.rowData <- read.counts.raw.rowData[chr.order, ]

  read.counts.raw <- read.counts.raw[chr.order, ]

  read.counts.raw.rowData$chr <- factor(read.counts.raw.rowData$chr, levels = unique(read.counts.raw.rowData$chr))

  # basic metrics
  read.counts.raw.colData$total.reads <- colSums(read.counts.raw) # total reads per cell
  read.counts.raw.rowData$total.reads <- rowSums(read.counts.raw) # total reads per probe
  read.counts.raw.rowData$median.reads <- apply(read.counts.raw, 1, median) # median reads per probe
  mean.reads <- round(mean(colMeans(read.counts.raw)), 2) # mean reads/cell/probe(or amplicon)

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
    colData = read.counts.raw.colData,
    rowData = read.counts.raw.rowData
  )

  tapestri.object <- .TapestriExperiment(sce)

  SingleCellExperiment::mainExpName(tapestri.object) <- "CNV"

  # save experiment metadata
  S4Vectors::metadata(tapestri.object)$sample.name <- "test object"
  S4Vectors::metadata(tapestri.object)$pipeline.panel.name <- "test object"
  S4Vectors::metadata(tapestri.object)$pipeline.version <- "test object"
  S4Vectors::metadata(tapestri.object)$number.of.cells <- 300
  S4Vectors::metadata(tapestri.object)$number.of.probes <- 304
  S4Vectors::metadata(tapestri.object)$date.h5.created <- "test object"
  S4Vectors::metadata(tapestri.object)$mean.reads.per.cell.per.probe <- as.character(mean.reads)

  # apply panel ID probe shortcuts
  tapestri.object@barcodeProbe <- barcodeProbe
  tapestri.object@grnaProbe <- grnaProbe

  tapestri.object <- getCytobands(tapestri.object)

  # move non-genomic probes to altExp slots
  tapestri.object <- moveNonGenomeProbes(tapestri.object)

  return(tapestri.object)
}

#' Create Example `TapestriExperiment`
#'
#' Creates a `TapestriExperiment` object for demonstration purposes,
#' which includes 240 probes across the genome, and 300 cells of 3 types.
#' Raw counts are generated randomly.
#' Type 1 has 75 cells, all XY, all diploid.
#' Type 2 has 100 cells, all XX, with 3 copies of chr 7, otherwise diploid.
#' Type 3 has 125 cells, all XY, with 1 copy of chr 1p, otherwise diploid.
#'
#' @importFrom stats rnorm
#'
#' @concept misc
#'
#' @return `TapestriExperiment` object.
#' @export
#'
#' @examples
#' tapExperiment <- newTapestriExperimentExample()
newTapestriExperimentExample <- function() {
  # panel ID
  barcodeProbe <- "dummyBCprobe"
  grnaProbe <- "dummyGRNAprobe"

  # probe metadata
  read.counts.raw.rowData <- S4Vectors::DataFrame(
    probe.id = paste0("probe_", 1:240),
    chr = rep(c(1:22, "X", "Y"), each = 10),
    arm = paste0("chr", gtools::mixedsort(c(rep(paste0(c(1:22, "X", "Y"), "p"), each = 5), rep(paste0(c(1:22, "X", "Y"), "q"), each = 5)))),
    start.pos = rep(1, 240),
    end.pos = rep(1, 240),
    row.names = paste0("probe_", 1:240)
  )

  read.counts.raw.rowData$chr <- factor(read.counts.raw.rowData$chr, levels = unique(read.counts.raw.rowData$chr))
  read.counts.raw.rowData$arm <- factor(read.counts.raw.rowData$arm, levels = unique(read.counts.raw.rowData$arm))

  n.probes <- 240

  # cell metadata
  read.counts.raw.colData <- S4Vectors::DataFrame(
    cell.barcode = as.character(paste0("cell_", 1:300)),
    test.cluster = c(rep("cellline1", 75), rep("cellline2", 100), rep("cellline3", 125)),
    row.names = as.character(paste0("cell_", 1:300))
  )

  # raw counts
  results.list <- list()
  for (i in 1:n.probes) {
    results.list[[i]] <- round(stats::rnorm(300, 100, 5), 0)
    names(results.list)[i] <- paste0("probe_", i)
  }

  read.counts.raw <- t(as.matrix(as.data.frame(results.list, row.names = NULL)))

  # Y counts
  read.counts.raw[which(read.counts.raw.rowData$chr == "Y"), c(1:75, 176:300)] <- 0

  # chr7 gain
  results.list <- list()
  for (i in 1:10) {
    results.list[[i]] <- round(stats::rnorm(100, 200, 5), 0)
  }

  read.counts.raw[which(read.counts.raw.rowData$chr == "7"), c(76:175)] <- t(as.matrix(as.data.frame(results.list, row.names = NULL)))

  # chr 1p loss
  results.list <- list()
  for (i in 1:5) {
    results.list[[i]] <- round(stats::rnorm(125, 50, 5), 0)
  }

  read.counts.raw[which(read.counts.raw.rowData$arm == "chr1p"), c(176:300)] <- t(as.matrix(as.data.frame(results.list, row.names = NULL)))


  # lentiviral probes
  ## skip for now

  # basic metrics
  read.counts.raw.colData$total.reads <- colSums(read.counts.raw) # total reads per cell
  read.counts.raw.rowData$total.reads <- rowSums(read.counts.raw) # total reads per probe
  read.counts.raw.rowData$median.reads <- apply(read.counts.raw, 1, median) # median reads per probe
  mean.reads <- round(mean(colMeans(read.counts.raw)), 2) # mean reads/cell/probe(or amplicon)

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = read.counts.raw),
    colData = read.counts.raw.colData,
    rowData = read.counts.raw.rowData
  )

  # construction
  tapestri.object <- .TapestriExperiment(sce)

  SingleCellExperiment::mainExpName(tapestri.object) <- "CNV"

  # save experiment metadata
  S4Vectors::metadata(tapestri.object)$sample.name <- "test object"
  S4Vectors::metadata(tapestri.object)$pipeline.panel.name <- "test object"
  S4Vectors::metadata(tapestri.object)$pipeline.version <- "test object"
  S4Vectors::metadata(tapestri.object)$number.of.cells <- 300
  S4Vectors::metadata(tapestri.object)$number.of.probes <- 240
  S4Vectors::metadata(tapestri.object)$date.h5.created <- "test object"
  S4Vectors::metadata(tapestri.object)$mean.reads.per.cell.per.probe <- as.character(mean.reads)

  # apply panel ID probe shortcuts
  tapestri.object@barcodeProbe <- barcodeProbe
  tapestri.object@grnaProbe <- grnaProbe

  # move non-genomic probes to altExp slots
  tapestri.object <- moveNonGenomeProbes(tapestri.object)


  # allele frequency
  af.matrix <- matrix(0, nrow = 75, ncol = 300)
  af.matrix[1:25, 1:75] <- 100
  af.matrix[26:50, 76:175] <- 100
  af.matrix[51:75, 175:300] <- 100

  colnames(af.matrix) <- colnames(tapestri.object)

  variant.metadata <- data.frame(variant.id = paste0("var_", 1:75))

  allele.frequency <- SingleCellExperiment::SingleCellExperiment(list(alleleFrequency = af.matrix),
    rowData = S4Vectors::DataFrame(variant.metadata),
    colData = S4Vectors::DataFrame(cell.barcode = colnames(tapestri.object))
  )

  allele.frequency <- .TapestriExperiment(allele.frequency)
  allele.frequency@barcodeProbe <- barcodeProbe
  allele.frequency@grnaProbe <- grnaProbe

  SingleCellExperiment::altExp(tapestri.object, "alleleFrequency", withDimnames = TRUE) <- allele.frequency

  return(tapestri.object)
}