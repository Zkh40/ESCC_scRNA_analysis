# Malignant epithelial t-SNE, stemness/EMT, CytoTRACE, ITH and Hallmark GSVA.
# Run from this folder: Rscript --vanilla 08_epithelial_state.R
source("00_paths.R")

align_lesion_metadata <- function (x, metadata) 
{
    stopifnot(!anyDuplicated(metadata$cell_id) && setequal(colnames(x), metadata$cell_id))
    m <- metadata[match(colnames(x), metadata$cell_id), , drop = FALSE]
    stopifnot(!anyNA(m$sample_id) && all(nzchar(m$sample_id)))
    m
}

lesion_mean_expression <- function (normalized_data, metadata) 
{
    m <- align_lesion_metadata(normalized_data, metadata)
    ids <- unique(m$sample_id)
    out <- vapply(ids, function(id) Matrix::rowMeans(normalized_data[, m$sample_id == id, drop = FALSE]), 
        numeric(nrow(normalized_data)))
    dimnames(out) <- list(rownames(normalized_data), ids)
    out
}

ith_by_lesion <- function (normalized_data, hvg_genes, metadata, config) 
{
    stopifnot(length(hvg_genes) == config$ith_hvg_n && !anyDuplicated(hvg_genes) && all(hvg_genes %in% 
        rownames(normalized_data)))
    m <- align_lesion_metadata(normalized_data, metadata)
    do.call(rbind, lapply(unique(m$sample_id), function(sid) {
        x <- as.matrix(normalized_data[hvg_genes, m$sample_id == sid, drop = FALSE])
        if (ncol(x) < 2) 
            return(data.frame(sample_id = sid, cells = ncol(x), pairs = 0L, undefined_pairs = 0L, ITH = NA_real_))
        r <- stats::cor(x, method = "pearson", use = "pairwise.complete.obs")
        v <- r[upper.tri(r, diag = FALSE)]
        good <- is.finite(v)
        stopifnot(all(good))
        data.frame(sample_id = sid, cells = ncol(x), pairs = length(v), undefined_pairs = sum(!good), 
            ITH = if (any(good)) 
                stats::IQR(v[good], type = 7)
            else NA_real_)
    }))
}

# ---- Analysis ----
library(Seurat)
# Frozen malignant-like selection with RNA PCA, counts and data. Selection is reviewed, not aneuploid alone.
obj <- readRDS(input("epithelial_malignant.rds"))
stopifnot(all(obj$copykat_malignancy_call == "MALIGNANT_LIKE"))
set.seed(1234)
fit <- Rtsne::Rtsne(Embeddings(obj, "pca")[, 1:25], dims = 2, pca = FALSE, perplexity = 30,
                    theta = 0.5, max_iter = 1000, check_duplicates = TRUE, partial_pca = FALSE,
                    is_distance = FALSE, Y_init = NULL, pca_center = TRUE, pca_scale = FALSE,
                    normalize = TRUE, stop_lying_iter = 250L, mom_switch_iter = 250L,
                    momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 1)
write_tab(data.frame(cell_id = colnames(obj), tSNE_1 = fit$Y[, 1], tSNE_2 = fit$Y[, 2]), "epithelial_tsne.tsv")
genes <- read_tab("resources/gene_sets.tsv")
for (nm in c("HCC10", "EMT_Mes29")) {
  features <- genes$gene[genes$signature == nm]
  obj <- AddModuleScore(obj, features = list(features), assay = "RNA", name = paste0(nm, "_"),
                        nbin = 24, ctrl = 100, seed = 1, search = FALSE)
}
meta <- data.frame(cell_id = colnames(obj), obj[[]])
write_tab(meta, "epithelial_scores.tsv")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- counts[Matrix::rowSums(counts > 0) >= 5, , drop = FALSE]
ct <- CytoTRACE::CytoTRACE(as.matrix(counts), ncores = 1) # Historical seed not recorded.
write_tab(data.frame(cell_id = names(ct$CytoTRACE), CytoTRACE = unname(ct$CytoTRACE)), "epithelial_CytoTRACE.tsv")
expr <- GetAssayData(obj, assay = "RNA", layer = "data")
hvg <- read_tab("resources/ith_malignant_hvg3000.tsv")$gene
write_tab(ith_by_lesion(expr, hvg, meta, list(ith_hvg_n = 3000)), "epithelial_ITH.tsv")
hallmark <- readRDS(input("hallmark50.rds")) # named list of exactly 50 frozen gene sets
stopifnot(length(hallmark) == 50)
means <- lesion_mean_expression(expr, meta) # mean of log RNA/data, unlike script 12
raw <- GSVA::gsva(GSVA::gsvaParam(means, hallmark, kcdf = "Gaussian", tau = 1,
                                maxDiff = TRUE, absRanking = FALSE), verbose = FALSE)
ann <- unique(meta[, c("sample_id", "IMMO")]); group <- ann$IMMO[match(colnames(raw), ann$sample_id)]
stopifnot(sum(group == "IM") == 12, sum(group == "MO") == 5, nrow(raw) == 50)
stats <- data.frame(pathway = rownames(raw), IM_median = apply(raw[, group == "IM"], 1, median),
                    MO_median = apply(raw[, group == "MO"], 1, median),
                    P = apply(raw, 1, function(x) wilcox.test(x[group == "IM"], x[group == "MO"], exact = FALSE)$p.value))
stats$BH_FDR <- p.adjust(stats$P, "BH")
write_tab(stats, "epithelial_Hallmark_statistics.tsv")
write_tab(data.frame(pathway = rownames(raw), raw), "epithelial_Hallmark_raw.tsv")
write_tab(data.frame(pathway = rownames(raw), t(scale(t(raw)))), "epithelial_Hallmark_row_z.tsv")
