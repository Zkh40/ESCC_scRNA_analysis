# Per-sample CopyKAT from original UMI with same-sample T/NK references.
# Run from this folder: Rscript --vanilla 07_epithelial_copykat.R
source("00_paths.R")


# ---- Analysis ----
matrices <- read_tab(input("matrix_manifest.tsv"))
cells <- read_tab(input("copykat_cells.tsv")) # sample_id, barcode_raw, role (query/reference)
for (sid in unique(cells$sample_id)) {
  umi <- Seurat::Read10X(matrices$matrix_dir[match(sid, matrices$sample_id)], gene.column = 2, unique.features = TRUE)
  m <- cells[cells$sample_id == sid, ]
  query <- m$barcode_raw[m$role == "query"]
  reference <- m$barcode_raw[m$role == "reference"]
  stopifnot(length(query) > 0, length(reference) > 0, !anyDuplicated(c(query, reference)),
            all(c(query, reference) %in% colnames(umi)))
  # CopyKAT 1.1.0 needs a base double matrix; preserve original UMI exactly.
  rawmat <- as.matrix(umi[, c(query, reference), drop = FALSE])
  storage.mode(rawmat) <- "double"
  native_dir <- normalizePath(output("."), winslash = "/")
  native_dir <- file.path(native_dir, paste0("copykat_", sid))
  stopifnot(!dir.exists(native_dir))
  dir.create(native_dir)
  original_wd <- getwd()
  setwd(native_dir)
  set.seed(42)
  ans <- copykat::copykat(rawmat = rawmat, id.type = "S", cell.line = "no", ngene.chr = 5,
                         min.gene.per.cell = 200, LOW.DR = 0.05, UP.DR = 0.1, win.size = 25,
                         KS.cut = 0.1, sam.name = sid, distance = "euclidean", norm.cell.names = reference,
                         output.seg = FALSE, plot.genes = TRUE, genome = "hg20", n.cores = 4)
  setwd(original_wd)
  pred <- ans$prediction
  pred <- pred[match(query, pred$cell.names), , drop = FALSE]
  pred$cell_id <- paste0(sid, "_", pred$cell.names)
  # not.defined remains unresolved; aneuploid/diploid are evidence for reviewed malignant/diploid-like labels.
  write_tab(pred, paste0(sid, "_copykat_query_predictions.tsv"))
}
