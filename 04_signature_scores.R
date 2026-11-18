# Score frozen gene panels on their original parent populations.
# Run from this folder: Rscript --vanilla 04_signature_scores.R
source("00_paths.R")

compare_scores <- function (df, comparisons, config) 
{
    stopifnot(!anyDuplicated(df$barcode) && !anyNA(df$IMMO) && all(df$IMMO %in% c("IM", "MO")))
    stopifnot(!anyDuplicated(comparisons$id))
    out <- lapply(seq_len(nrow(comparisons)), function(i) {
        z <- comparisons[i, ]
        d <- df
        if ("cluster" %in% names(z)) {
            d <- df[df$cluster == z$cluster, , drop = FALSE]
        }
        x <- d[[z$score]][d$IMMO == "IM"]
        y <- d[[z$score]][d$IMMO == "MO"]
        stopifnot(length(x) > 0 && length(y) > 0 && all(is.finite(c(x, y))))
        qs <- function(v) unname(stats::quantile(v, c(0.25, 0.5, 0.75)))
        qx <- qs(x)
        qy <- qs(y)
        data.frame(id = z$id, family = z$family, n_IM = length(x), n_MO = length(y), Q1_IM = qx[1], median_IM = qx[2], 
            Q3_IM = qx[3], Q1_MO = qy[1], median_MO = qy[2], Q3_MO = qy[3], cliff_delta = 2 * (sum(rank(c(x, 
                y), ties.method = "average")[seq_along(x)]) - length(x) * (length(x) + 1)/2)/(length(x) * 
                length(y)) - 1, P = stats::wilcox.test(x, y, alternative = "two.sided", exact = FALSE, 
                correct = TRUE)$p.value)
    })
    out <- do.call(rbind, out)
    out$FDR <- NA_real_
    for (f in unique(out$family)) {
        if (is.na(f) || f == "NOT_PROVIDED") 
            next
        expected <- config$families[[f]]
        stopifnot(length(expected) > 0 && !anyDuplicated(expected) && setequal(out$id[out$family == f], 
            expected))
        j <- which(out$family == f)
        out$FDR[j] <- stats::p.adjust(out$P[j], "BH")
    }
    out
}

# ---- Analysis ----
library(Seurat)
genes <- read_tab("resources/gene_sets.tsv")
panels <- list(TNK = c("Chu14", "Ttex", "Tpex", "TRM", "Treg"),
               Mac = c("Phagosome", "M1_like", "IFNG_response", "Antigen_presentation", "Angiogenesis", "MHC_I", "MHC_II"),
               DC = c("DC_activation", "DC_tolerance"))
# These supplied objects preserve the historical scoring parents: clean TNK; Mac 14365; DC C13 753.
# A final cluster subset must not be scored as if it were its original parent population.
for (population in names(panels)) {
  obj <- readRDS(input(paste0(population, "_score_parent.rds")))
  for (nm in panels[[population]]) {
    features <- genes$gene[genes$signature == nm]
    stopifnot(length(features) > 0, all(features %in% rownames(obj[["RNA"]])))
    obj <- AddModuleScore(obj, features = list(features), assay = "RNA", slot = "data",
                          name = paste0(nm, "_"), nbin = 24, ctrl = 100, seed = 1, search = FALSE)
  }
  meta <- data.frame(barcode = colnames(obj), obj[[]])
  # AddModuleScore appends 1; use the signature names in the comparison table.
  for (nm in panels[[population]]) meta[[nm]] <- meta[[paste0(nm, "_1")]]
  write_tab(meta, paste0(population, "_scores.tsv"))
  if (population == "TNK") {
    comparisons <- read_tab("resources/tnk_comparisons_60.tsv")
    config <- list(families = list(Ttex_Tpex_TRM_60 = comparisons$id))
    write_tab(compare_scores(meta, comparisons, config), "TNK_60_score_comparisons.tsv")
  }
}
# Historical cell-level tests are exploratory. Treg has no reconstructed FDR.
# MHC FDR requires all four original KEGG/GO comparisons; do not BH-adjust the two GO panels alone.
