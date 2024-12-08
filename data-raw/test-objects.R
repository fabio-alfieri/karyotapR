tap.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.dna.h5"
h5.filename <- tap.file

bam.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.tube1.cells.bam"

barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

#full set
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3"),
                             sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG", "ACTCTTGCTGTGGCATTTTC"))

#Troubleshooting object, gNC is random sequence
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
                             sequences = c("TGGATATATGGACCGCATTG", "CGACTGTTCCCAAATTGTAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG"))

#Troubleshooting object, all bcs are random sequences
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2"),
           sequences = c("GATGTAACCATATACTTACG", "CTGGATCTTCTCCCGCGAAT", "TTTAACCCTCACCAACTACG", "AGATTTGAGGTAAACCAAAT"))

parsed.bc <- countBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "GRNA")

countBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "barcode")


bam.matches <- bam.filter$tag[sequence.matches[[x]]]

coldata.labels <- c("g7x1", "gNC")

library(devtools)
load_all()
tap.file <- "~/Documents/data-and-analyses/tapestri/data/20220331-CO293/20220331-CO293-reference20220419.dna.h5"

tap <- createTapestriExperiment(tap.file, "CO293")
tap <- countBarcodedReads(tap, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "grna")
colData(tap)$sample.grna <- callSampleLables(tap, coldata.labels = c(coldata.labels, "g9p21L2", "g9p21R2"))
tap1 <- runPCA(tap, sd.min.threshold = 35)
PCAKneePlot(tap1)
tap2 <- runUMAP(tap1, input.dims = 1:3)
reducedDimPlot(tap2, "umap", group.label = "sample.grna")
umap.dbscan <- dbscan::dbscan(reducedDim(altExp(tap2, "alleleFrequency"), "UMAP"), eps = 0.8)
colData(tap2)$cluster <- as.factor(umap.dbscan$cluster)
reducedDimPlot(tap2, "umap", group.label = "cluster")
reducedDimPlot(tap2, "umap", group.label = "sample.grna")


library(tidyverse)
load_all()
tap.file <- "~/Documents/data-and-analyses/tapestri/data/Exp003-20220808-CO293/run20220808-02-panelCO293-ref20220419.dna.h5"
bam.file <- "~/Documents/data-and-analyses/tapestri/data/Exp003-20220808-CO293/run20220808-02-panelCO293-ref20220419.tube1.cells.bam"
barcode.lookup <- data.frame(ids = c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3"), sequences = c("TGGATATATGGACCGCATTG", "ACGGAGGCTAAGCGTCGCAA", "TATTTACAGGGACAATACCG", "CGGTAGAATAAGCTGTACCG", "ACTCTTGCTGTGGCATTTTC"))
label.features <- c("g7x1", "gNC", "g9p21L2", "g9p21R2", "g7x3")

exp3 <- createTapestriExperiment(tap.file, "CO293")
exp3 <- countBarcodedReads(exp3, bam.file = bam.file, barcode.lookup = barcode.lookup, probe = "grna")
exp3 <- callSampleLables(exp3, input.features = c("g7x3", "gNC"), output.feature = "sample.grna")
exp3 <- runPCA(exp3, sd.min.threshold = 40)
PCAKneePlot(exp3)
exp3 <- runUMAP(exp3, pca.dims = 1:2)
reducedDimPlot(exp3, dim.reduction = "umap")
exp3 <- runClustering(exp3, eps = 0.9)
reducedDimPlot(exp3, dim.reduction = "umap", group.label = "cluster")
reducedDimPlot(exp3, dim.reduction = "umap", group.label = "sample.grna")
forcats::fct_count(colData(exp3)$cluster)
exp3.subset <- exp3[,colData(exp3)$cluster %in% 1:3]
reducedDimPlot(exp3.subset, dim.reduction = "umap", group.label = "cluster")
colData(exp3.subset)$cluster <- forcats::fct_recode(colData(exp3.subset)$cluster, RPE1 = "1", cellline1 = "2", cellline2 = "3")
colData(exp3.subset)$cluster <- forcats::fct_drop(colData(exp3.subset)$cluster)
forcats::fct_count(colData(exp3.subset)$cluster)
summary(rowData(exp3.subset)$median.reads > 0)
exp3.subset <- exp3.subset[rowData(exp3.subset)$median.reads > 0,] #1 probe filtered for exp 3
exp3.subset <- calcNormCounts(exp3.subset)

control.copy.number <- generateControlCopyNumberTemplate(exp3.subset, sample.feature.label = "RPE1")
control.copy.number["chr10q", "copy.number"] <- 3
exp3.subset <- calcCopyNumber(exp3.subset, control.copy.number = control.copy.number, sample.feature = "cluster")

exp3.subset <- calcSmoothCopyNumber(exp3.subset, method = "median")
exp3.subset <- calcGMMCopyNumber(exp3.subset, cell.barcodes = colnames(exp3.subset), control.copy.number = control.copy.number, model.components = 1:4)

exp3.subset <- calcSmoothCopyNumber(exp3.subset, method = "weighted.median", control.copy.number = control.copy.number)
exp3.subset <- calcGMMCopyNumber(exp3.subset, cell.barcodes = colnames(exp3.subset), control.copy.number = control.copy.number, model.components = 1:4, 
                                 smoothing.method = "weighted.median")

assayHeatmap(exp3.subset, assay = "copyNumber", split.col.by = "chr", split.row.by = "cluster", annotate.row.by = "sample.grna", color.preset = "copy.number.denoise")
assayHeatmap(exp3.subset, alt.exp = "smoothedCopyNumberByChr", assay = "discreteCopyNumber", split.row.by = "cluster", annotate.row.by = "sample.grna", color.preset = "copy.number")
assayBoxPlot(exp3.subset, alt.exp = "chrYCounts", split.features = TRUE, split.x.by = "cluster")

getTidyData(exp3.subset)
getTidyData(exp3.subset, alt.exp = "alleleFrequency")
getTidyData(exp3.subset, assay = "copyNumber")
getTidyData(exp3.subset, alt.exp = "chrYCounts")
getTidyData(exp3.subset, alt.exp = "smoothedCopyNumberByChr")


# readme test

example.exp <- newTapestriExperimentExample()
assayHeatmap(example.exp, split.col.by = "arm", split.row.by = "test.cluster", annotate.row.by = "test.cluster")

example.exp <- runPCA(example.exp)
PCAKneePlot(example.exp)
example.exp <- runUMAP(example.exp, pca.dims = 1:2)
reducedDimPlot(example.exp, dim.reduction = "umap")
example.exp <- runClustering(example.exp, eps = 0.9)
reducedDimPlot(example.exp, dim.reduction = "umap", group.label = "cluster")
colData(example.exp)$cluster <- forcats::fct_recode(colData(example.exp)$cluster, cellline1 = "1", cellline2 = "2", cellline3 = "3")

example.exp <- calcNormCounts(example.exp)
control.copy.number <- generateControlCopyNumberTemplate(example.exp, sample.feature.label = "cellline3", copy.number = 2)
example.exp <- calcCopyNumber(example.exp, control.copy.number = control.copy.number, sample.feature = "cluster")
example.exp <- calcSmoothCopyNumber(example.exp)


assayHeatmap(example.exp, assay = "copyNumber", split.col.by = "arm", split.row.by = "test.cluster", annotate.row.by = "test.cluster", color.preset = "copy.number")
assayHeatmap(example.exp, alt.exp = "smoothedCopyNumberByArm", assay = "discreteCopyNumber", split.row.by = "test.cluster", annotate.row.by = "test.cluster", color.preset = "copy.number")
assayBoxPlot(example.exp, alt.exp = "chrYCounts", split.features = TRUE, split.x.by = "test.cluster")



