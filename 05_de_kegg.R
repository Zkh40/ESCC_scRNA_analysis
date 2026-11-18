# Cell-level IM versus MO DE, T/DC ORA and macrophage ranked GSEA.
# Run from this folder: Rscript --vanilla 05_de_kegg.R
source("00_paths.R")

map_entrez_ranks <- function (deg, mapping) 
{
    stopifnot(!anyDuplicated(deg$SYMBOL) && all(is.finite(deg$avg_log2FC)))
    a <- unique(mapping[, c("SYMBOL", "ENTREZID")])
    a$ENTREZID <- as.character(a$ENTREZID)
    a <- a[a$SYMBOL %in% deg$SYMBOL, , drop = FALSE]
    n <- table(a$SYMBOL[!is.na(a$ENTREZID) & nzchar(a$ENTREZID)])
    a$stat <- deg$avg_log2FC[match(a$SYMBOL, deg$SYMBOL)]
    a$status <- ifelse(is.na(a$ENTREZID) | !nzchar(a$ENTREZID), "UNMAPPED", ifelse(n[a$SYMBOL] > 1, "DROP_AMBIGUOUS_SYMBOL", 
        "CANDIDATE"))
    u <- a[a$status == "CANDIDATE", , drop = FALSE]
    u <- u[order(-abs(u$stat), u$SYMBOL, u$ENTREZID), , drop = FALSE]
    kept <- u[!duplicated(u$ENTREZID), , drop = FALSE]
    keys <- paste(kept$SYMBOL, kept$ENTREZID)
    j <- which(a$status == "CANDIDATE")
    a$status[j] <- ifelse(paste(a$SYMBOL[j], a$ENTREZID[j]) %in% keys, "KEPT", "DROP_DUPLICATE_ENTREZ")
    kept <- kept[order(-kept$stat, kept$SYMBOL, kept$ENTREZID), , drop = FALSE]
    stats::setNames(kept$stat, kept$ENTREZID)
}

# ---- Analysis ----
library(Seurat)
library(clusterProfiler)
term2gene <- read_tab(input("TERM2GENE_hsa.tsv"))
term2name <- read_tab(input("TERM2NAME_hsa.tsv"))
jobs <- read_tab(input("de_jobs.tsv")) # id, profile (T/DC/Mac), object, feature_symbols
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  obj <- readRDS(job$object) # exact selected cells; metadata IMMO
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")
  im <- colnames(obj)[obj$IMMO == "IM"]; mo <- colnames(obj)[obj$IMMO == "MO"]
  min_pct <- if (job$profile == "DC") 0.25 else 0.1
  lfc <- if (job$profile == "DC") 0.25 else 0
  features <- rownames(expr)
  if (job$profile == "Mac") {
    pim <- Matrix::rowMeans(expr[, im, drop = FALSE] > 0)
    pmo <- Matrix::rowMeans(expr[, mo, drop = FALSE] > 0)
    features <- features[pmax(pim, pmo) >= 0.1]
  }
  assay <- SeuratObject::CreateAssayObject(data = expr)
  args <- list(object = assay, cells.1 = im, cells.2 = mo, features = features,
               slot = "data", fc.slot = "data", test.use = "wilcox", min.pct = min_pct,
               logfc.threshold = lfc, min.diff.pct = -Inf, only.pos = FALSE,
               max.cells.per.ident = Inf, pseudocount.use = 1, base = 2,
               fc.name = "avg_log2FC", densify = FALSE, verbose = FALSE)
  if (job$profile == "Mac") { set.seed(20261031); args$random.seed <- 20261031 }
  deg <- do.call(Seurat::FindMarkers, args)
  deg$SYMBOL <- rownames(deg)
  deg$BH_FDR <- p.adjust(deg$p_val, "BH")
  write_tab(deg, paste0(job$id, "_DE.tsv"))
  if (job$profile %in% c("T", "DC")) {
    # All original RNA features, not only DE rows. DC uses full Myeloid features.
    symbols <- read_tab(job$feature_symbols)$SYMBOL
    ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = symbols,
                                column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    universe <- intersect(unique(na.omit(as.character(ids))), unique(as.character(term2gene[[2]])))
    for (group in c("IM", "MO")) {
      selected <- deg$SYMBOL[deg$BH_FDR < 0.05 &
                   if (group == "IM") deg$avg_log2FC >= 0.25 else deg$avg_log2FC <= -0.25]
      selected_ids <- unique(na.omit(as.character(ids[selected])))
      ora <- enricher(selected_ids, universe = universe, TERM2GENE = term2gene, TERM2NAME = term2name,
                      minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
      write_tab(as.data.frame(ora), paste0(job$id, "_", group, "_KEGG_ORA.tsv"))
    }
  } else {
    mapping <- read_tab(input(paste0(job$id, "_symbol_entrez.tsv"))) # exact frozen mapping, SYMBOL/ENTREZID
    ranks <- map_entrez_ranks(deg, mapping)
    set.seed(20261031)
    fit <- GSEA(geneList = ranks, TERM2GENE = unique(term2gene), TERM2NAME = unique(term2name),
                minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH", pvalueCutoff = 1,
                eps = 0, seed = FALSE, by = "fgsea", nPermSimple = 10000,
                BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    write_tab(as.data.frame(fit), paste0(job$id, "_KEGG_GSEA.tsv"))
  }
}
# These reproduce exploratory cell-level tests; they are not patient-level inference.
