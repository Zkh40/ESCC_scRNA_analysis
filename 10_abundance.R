# Lesion-level state fractions and Spearman correlations; denominators differ by question.
# Run from this folder: Rscript --vanilla 10_abundance.R
source("00_paths.R")


# ---- Analysis ----
metadata <- read_tab("resources/samples.tsv")
global_states <- read_tab("resources/abundance_states.tsv")$state
for (scope in c("global", "targeted")) {
  tab <- read_tab(input(paste0(scope, "_abundance_counts.tsv"))) # sample_id, state, count, denominator
  states <- if (scope == "global") global_states else c("Mac_C1", "TNK_C2")
  # Global: TNK states / clean TNK INCLUDING NK; myeloid states / all myeloid.
  # Targeted: TNK_C2 / clean T EXCLUDING NK; Mac_C1 / all myeloid.
  # Long tables must include all measured zeros. Never substitute total captured cells.
  stopifnot(nrow(tab) == nrow(metadata) * length(states),
            !anyDuplicated(paste(tab$sample_id, tab$state)),
            all(tab$denominator > 0), all(tab$count >= 0), all(tab$count <= tab$denominator))
  x <- matrix(NA_real_, nrow(metadata), length(states), dimnames = list(metadata$sample_id, states))
  x[cbind(match(tab$sample_id, metadata$sample_id), match(tab$state, states))] <- tab$count / tab$denominator
  stopifnot(!anyNA(x), all(apply(x, 2, sd) > 0))
  pairs <- combn(states, 2)
  result <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
    a <- pairs[1, i]; b <- pairs[2, i]
    z <- cor.test(x[, a], x[, b], method = "spearman", exact = FALSE)
    data.frame(state_A = a, state_B = b, rho = unname(z$estimate), P = z$p.value)
  }))
  result$FDR <- if (scope == "global") p.adjust(result$P, "BH") else NA_real_
  write_tab(result, paste0(scope, "_abundance_correlations.tsv"))
  write_tab(data.frame(sample_id = rownames(x), x), paste0(scope, "_abundance_fractions.tsv"))
}
# Global BH uses all 703 pairs across all 17 lesions (12 IM / 5 MO).
# IM/MO-separated correlations are descriptive; repeated lesions are not independent patients.
