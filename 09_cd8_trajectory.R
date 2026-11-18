# CD8 Monocle2: new-tree estimation OR reroot the historical tree; separate processes.
# Run from this folder: Rscript --vanilla 09_cd8_trajectory.R
source("00_paths.R")

trajectory_stage1 <- function (counts, cluster_labels, config) 
{
    stopifnot(as.character(utils::packageVersion("monocle")) == "2.34.0")
    scope <- c("C2", "C3", "C5", "C8", "C11", "C16", "C19")
    stopifnot(identical(colnames(counts), as.character(cluster_labels$barcode)) && !anyDuplicated(cluster_labels$barcode))
    stopifnot(setequal(cluster_labels$cluster, scope))
    stopifnot(!anyDuplicated(rownames(counts)) && !is.null(rownames(counts)))
    if (!is.null(config$seed)) 
        set.seed(config$seed)
    pd <- data.frame(cluster = factor(cluster_labels$cluster, levels = scope), row.names = cluster_labels$barcode)
    fd <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
    cds <- monocle::newCellDataSet(counts, phenoData = Biobase::AnnotatedDataFrame(pd), featureData = Biobase::AnnotatedDataFrame(fd), 
        lowerDetectionLimit = 0.1, expressionFamily = VGAM::negbinomial.size())
    cds <- BiocGenerics::estimateSizeFactors(cds)
    cds <- BiocGenerics::estimateDispersions(cds)
    cds <- monocle::detectGenes(cds, min_expr = 0.1)
    expressed <- rownames(Biobase::fData(cds))[Biobase::fData(cds)$num_cells_expressed >= 10]
    diff <- monocle::differentialGeneTest(cds[expressed, ], fullModelFormulaStr = "~cluster", cores = config$cores)
    diff$gene <- rownames(diff)
    ranked <- diff[which(diff$qval < 0.01), , drop = FALSE]
    ranked <- ranked[order(ranked$qval, ranked$gene), , drop = FALSE]
    genes <- head(ranked$gene, 400)
    stopifnot(length(genes) == 400 && !anyDuplicated(genes))
    cds <- monocle::setOrderingFilter(cds, genes)
    cds <- monocle::reduceDimension(cds, reduction_method = "DDRTree")
    cds <- orderCells_compat(cds)
    list(cds = cds, ordering_genes = ranked, differential = diff, seed = config$seed)
}

# ---- Analysis ----
# Use the supplied Monocle2-specific library profile (see software_versions.tsv).
source("functions/monocle2_compat.R")
mode <- Sys.getenv("ESCC_TRAJECTORY_MODE", unset = "reroot")
if (mode == "estimate") {
  counts <- readRDS(input("cd8_counts.rds"))
  labels <- read_tab(input("cd8_labels.tsv")) # barcode, cluster; columns in exact count order
  config <- list(seed = NULL, cores = 6)
  z <- trajectory_stage1(counts, labels, config)
  write_tab(z$differential, "CD8_ordering_gene_tests.tsv")
  saveRDS(z$cds, output("CD8_new_tree.rds"))
  # Historical DDRTree seed was not recorded. New State numbers/node IDs need fresh interpretation.
} else {
  stopifnot(mode == "reroot")
  cds <- readRDS(input("cd8_historical_tree.rds"))
  # State3 / Y_78 is defined only on this historical tree, never on CD8_new_tree.rds.
  projection <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  nodes <- colnames(monocle::reducedDimK(cds))[as.integer(projection)]
  pd <- Biobase::pData(cds)
  stopifnot(identical(rownames(projection), colnames(cds)),
            "Y_78" %in% nodes[pd$State == 3], igraph::degree(monocle::minSpanningTree(cds), v = "Y_78") == 1)
  cds <- fixed_order_cells(cds, root_state = 3)
  root_cell <- cds@auxOrderingData[["DDRTree"]]$root_cell
  stopifnot(nodes[match(root_cell, colnames(cds))] == "Y_78")
  write_tab(data.frame(barcode = rownames(Biobase::pData(cds)), Biobase::pData(cds)), "CD8_historical_pseudotime.tsv")
}
