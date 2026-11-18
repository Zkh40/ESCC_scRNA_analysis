# Global normalization, Harmony, clustering and broad annotation.
# Run from this folder: Rscript --vanilla 02_integrate_annotate.R
source("00_paths.R")


# ---- Analysis ----
library(Seurat)
obj <- readRDS(output("qc_singlets.rds"))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
cc <- Seurat::cc.genes.updated.2019
obj <- CellCycleScoring(obj, s.features = cc$s.genes, g2m.features = cc$g2m.genes, set.ident = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = NULL)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, seed.use = 12345)
set.seed(12345)
obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "pca", dims.use = 1:30)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.3, random.seed = 12345)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, min.dist = 0.2,
               n.neighbors = 30, metric = "cosine", seed.use = 42)
# Historical cluster names are meaningful only after marker-based review of a new clustering.
write_tab(data.frame(barcode = colnames(obj), obj[[]]), "global_clusters.tsv")
pdf(output("global_clusters.pdf")); print(DimPlot(obj, label = TRUE)); dev.off()
mapping <- read_tab(input("reviewed_broad_labels.tsv")) # cluster, label; see supplied historical reference
stopifnot(!anyDuplicated(mapping$cluster), setequal(as.character(unique(obj$seurat_clusters)), as.character(mapping$cluster)))
obj$broad_type <- mapping$label[match(as.character(obj$seurat_clusters), as.character(mapping$cluster))]
saveRDS(obj, output("atlas_annotated.rds"))
