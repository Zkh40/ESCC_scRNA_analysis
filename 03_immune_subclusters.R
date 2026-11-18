# Separate TNK and myeloid reclustering; review labels before downstream use.
# Run from this folder: Rscript --vanilla 03_immune_subclusters.R
source("00_paths.R")


# ---- Analysis ----
library(Seurat)
atlas <- readRDS(output("atlas_annotated.rds"))
for (lineage in c("TNK", "Myeloid")) {
  broad <- if (lineage == "TNK") "T/NK" else "Myeloid"
  obj <- subset(atlas, cells = colnames(atlas)[atlas$broad_type == broad])
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, seed.use = 42)
  obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca", dims.use = 1:20)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, k.param = 20)
  obj <- FindClusters(obj, resolution = 1.2, random.seed = 1234)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, n.neighbors = 30, min.dist = 0.3, seed.use = 1234)
  obj$cluster <- paste0("C", obj$seurat_clusters)
  write_tab(data.frame(barcode = colnames(obj), obj[[]]), paste0(lineage, "_clusters.tsv"))
  labels <- read_tab(input("reviewed_immune_labels.tsv")) # lineage, cluster, label, retained
  labels <- labels[labels$lineage == lineage, ]
  # Mark rejected TNK clusters explicitly in the reviewed table. Historical exclusions: C13/C15/C22/C23.
  stopifnot(!anyDuplicated(labels$cluster), all(obj$cluster %in% labels$cluster))
  obj$label <- labels$label[match(obj$cluster, labels$cluster)]
  retained <- labels$retained[match(obj$cluster, labels$cluster)]
  obj <- subset(obj, cells = colnames(obj)[retained])
  saveRDS(obj, output(paste0(lineage, "_clean.rds")))
}
