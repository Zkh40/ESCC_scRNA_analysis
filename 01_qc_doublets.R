# Original UMI → per-sample decontX → QC → sample-split doublet removal.
# Run from this folder: Rscript --vanilla 01_qc_doublets.R
source("00_paths.R")


# ---- Analysis ----
library(Seurat)
library(Matrix)
samples <- read_tab("resources/samples.tsv")
matrices <- read_tab(input("matrix_manifest.tsv")) # sample_id, matrix_dir
stopifnot(setequal(samples$sample_id, matrices$sample_id))
hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
corrected <- metadata <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  sid <- samples$sample_id[i]
  umi <- Read10X(matrices$matrix_dir[match(sid, matrices$sample_id)], gene.column = 2, unique.features = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = umi))
  sce <- celda::decontX(sce, seed = 12345)
  contamination <- SummarizedExperiment::colData(sce)$decontX_contamination
  keep <- contamination <= 0.2
  x <- SummarizedExperiment::assay(sce, "decontXcounts")[, keep, drop = FALSE]
  raw_barcodes <- colnames(x)
  colnames(x) <- paste0(sid, "_", raw_barcodes)
  corrected[[i]] <- as(x, "sparseMatrix")
  m <- samples[rep(i, sum(keep)), , drop = FALSE]
  rownames(m) <- colnames(x)
  m$barcode_raw <- raw_barcodes
  m$orig.ident <- sid
  m$decontX_contamination <- contamination[keep]
  metadata[[i]] <- m
}
stopifnot(all(vapply(corrected, function(x) identical(rownames(x), rownames(corrected[[1]])), logical(1))))
obj <- CreateSeuratObject(do.call(cbind, corrected), min.cells = 0, min.features = 0)
obj <- AddMetaData(obj, do.call(rbind, metadata))
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.HB"]] <- PercentageFeatureSet(obj, features = intersect(hb_genes, rownames(obj)))
obj <- subset(obj, subset = nFeature_RNA >= 500 & nFeature_RNA <= 7500 & percent.mt < 20 & percent.HB < 5)
set.seed(12345)
sce <- scDblFinder::scDblFinder(as.SingleCellExperiment(obj, assay = "RNA"),
                              samples = "sample_id", multiSampleMode = "split")
stopifnot(identical(colnames(sce), colnames(obj)))
obj$scDblFinder.class <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
obj <- subset(obj, subset = scDblFinder.class == "singlet")
saveRDS(obj, output("qc_singlets.rds"))
# Original matrices stay external and unchanged; CopyKAT uses them, not corrected counts.
