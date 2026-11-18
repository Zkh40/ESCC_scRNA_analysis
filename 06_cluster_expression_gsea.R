# Within-cluster expression enrichment; pooled IM/MO, not differential GSEA.
# Run from this folder: Rscript --vanilla 06_cluster_expression_gsea.R
source("00_paths.R")

cluster_expression_gsea <- function (expr, cell_ids, signature, pathway = "signature", seed = 20261107L) 
{
    x <- methods::as(expr[, cell_ids, drop = FALSE], "dgCMatrix")
    x@x <- expm1(x@x)
    means <- Matrix::rowMeans(x)
    keep <- means > 0
    means <- means[keep]
    names(means) <- rownames(x)[keep]
    stats <- means[order(-means, names(means), method = "radix")]
    overlap <- signature[signature %in% names(stats)]
    if (!length(overlap)) 
        stop("signature has no overlap with expressed genes")
    if (length(overlap) < 15L || length(overlap) > 500L) 
        stop("overlapping signature size outside prespecified 15..500")
    if (length(overlap) >= length(stats)) 
        stop("signature must be smaller than ranked universe")
    set.seed(as.integer(seed))
    res <- fgsea::fgseaMultilevel(pathways = setNames(list(signature), pathway), stats = stats, minSize = 15L, 
        maxSize = 500L, eps = 0, nPermSimple = 10000L, sampleSize = 101L, gseaParam = 1, scoreType = "std", 
        BPPARAM = BiocParallel::SerialParam(progressbar = FALSE))
    if (nrow(res) != 1L || !all(is.finite(c(res$pval, res$NES, res$ES))) || res$pval[1] < 0 || res$pval[1] > 
        1) 
        stop("fgsea did not return a valid result")
    lead <- res$leadingEdge[[1L]]
    rank_index <- match(lead, names(stats))
    leading_edge <- data.frame(gene = lead, rank = rank_index, mean_expression = unname(stats[rank_index]), 
        stringsAsFactors = FALSE)
    leading_edge <- leading_edge[order(leading_edge$rank), , drop = FALSE]
    membership <- data.frame(gene = signature, in_RNA = signature %in% rownames(expr), expressed = signature %in% 
        names(stats), leading_edge = signature %in% lead, stringsAsFactors = FALSE)
    result <- data.frame(pathway = res$pathway, size = res$size, ES = res$ES, NES = res$NES, rawP = res$pval, 
        log2err = res$log2err, stringsAsFactors = FALSE)
    result$precision_status <- ifelse(is.finite(result$log2err), "ESTIMATED", "REVIEW_FGSEA_WARNING")
    list(result = result, ranks = stats, leading_edge = leading_edge, membership = membership)
}

# ---- Analysis ----
obj <- readRDS(input("TNK_expression_parent.rds"))
expr <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
genes <- read_tab("resources/gene_sets.tsv")
jobs <- read_tab(input("expression_gsea_jobs.tsv")) # id, cluster, signature
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  cell_ids <- colnames(obj)[obj$cluster == job$cluster]
  signature <- genes$gene[genes$signature == job$signature]
  z <- cluster_expression_gsea(expr, cell_ids, signature, pathway = job$signature, seed = 20261107)
  write_tab(z$result, paste0(job$id, "_expression_GSEA.tsv"))
  write_tab(z$leading_edge, paste0(job$id, "_leading_edge.tsv"))
}
