# Within-patient lesion similarity, median per patient, exact 28-allocation comparison.
# Run from this folder: Rscript --vanilla 12_transcript_similarity.R
source("00_paths.R")

similarity_lesion_profiles <- function (corrected_counts, metadata, shared_hvg2000, config) 
{
    stopifnot(length(dim(corrected_counts)) == 2L && all(dim(corrected_counts) > 0))
    stopifnot(!is.null(rownames(corrected_counts)) && !is.null(colnames(corrected_counts)) && !anyDuplicated(rownames(corrected_counts)) && 
        !anyDuplicated(colnames(corrected_counts)))
    stopifnot(!anyDuplicated(metadata$cell_id) && setequal(metadata$cell_id, colnames(corrected_counts)))
    stopifnot(all(is.finite(corrected_counts)) && all(corrected_counts >= 0))
    stopifnot(length(shared_hvg2000) > 1L && !anyDuplicated(shared_hvg2000) && all(shared_hvg2000 %in% 
        rownames(corrected_counts)))
    stopifnot(length(shared_hvg2000) == 2000L)
    m <- metadata[match(colnames(corrected_counts), metadata$cell_id), , drop = FALSE]
    stopifnot(!anyNA(m[, c("sample_id", "patient_id", "IMMO")]) && all(nzchar(m$sample_id)))
    ann <- unique(m[, c("sample_id", "patient_id", "IMMO")])
    stopifnot(!anyDuplicated(ann$sample_id))
    totals <- if (inherits(corrected_counts, "Matrix")) 
        Matrix::colSums(corrected_counts)
    else colSums(corrected_counts)
    stopifnot(all(is.finite(totals)) && all(totals > 0))
    ids <- unique(m$sample_id)
    profiles <- matrix(NA_real_, length(shared_hvg2000), length(ids), dimnames = list(shared_hvg2000, 
        ids))
    for (sid in ids) {
        j <- which(m$sample_id == sid)
        if (length(j) < config$min_cells) 
            next
        x <- as.matrix(corrected_counts[shared_hvg2000, j, drop = FALSE])
        nonlog <- sweep(x, 2, config$scale_factor/totals[j], "*")
        profiles[, sid] <- log1p(rowMeans(nonlog))
    }
    list(profiles = profiles, lesions = ann, cell_counts = setNames(as.integer(table(factor(m$sample_id, 
        levels = ids))), ids))
}

similarity_within_patient_pairs <- function (profile_result) 
{
    ann <- profile_result$lesions
    p <- profile_result$profiles
    sizes <- profile_result$cell_counts
    rows <- list()
    for (pid in unique(ann$patient_id)) {
        a <- ann[ann$patient_id == pid, , drop = FALSE]
        stopifnot(length(unique(a$IMMO)) == 1L)
        stopifnot(nrow(a) >= 2L)
        pairs <- combn(a$sample_id, 2)
        for (i in seq_len(ncol(pairs))) {
            lhs <- pairs[1, i]
            rhs <- pairs[2, i]
            x <- p[, lhs]
            y <- p[, rhs]
            finite <- all(is.finite(x)) && all(is.finite(y))
            reason <- if (!finite) 
                "BELOW_CELL_THRESHOLD_OR_NONFINITE"
            else if (var(x) == 0 || var(y) == 0) 
                "ZERO_VARIANCE_VECTOR"
            else ""
            r <- if (nzchar(reason)) 
                NA_real_
            else cor(x, y, method = "pearson")
            rows[[length(rows) + 1L]] <- data.frame(patient_id = pid, IMMO = a$IMMO[1], lesion_A = lhs, 
                lesion_B = rhs, n_cells_A = unname(sizes[lhs]), n_cells_B = unname(sizes[rhs]), r = r, 
                eligibility_status = if (is.finite(r)) 
                  "ELIGIBLE"
                else "NOT_ESTIMABLE", exclusion_reason = reason)
        }
    }
    do.call(rbind, rows)
}

similarity_patient_summary <- function (pair_results, config) 
{
    do.call(rbind, lapply(unique(pair_results$patient_id), function(pid) {
        d <- pair_results[pair_results$patient_id == pid, , drop = FALSE]
        good <- is.finite(d$r)
        stopifnot(length(unique(d$IMMO)) == 1L)
        data.frame(patient_id = pid, IMMO = d$IMMO[1], patient_r = if (any(good)) 
            median(d$r[good])
        else NA_real_, eligible_pair_n = sum(good), expected_pair_n = nrow(d), aggregation = "median_of_eligible_pair_r")
    }))
}

similarity_exact_permutation <- function (patient_results, config) 
{
    stopifnot(!anyDuplicated(patient_results$patient_id) && nrow(patient_results) == 8L && sum(patient_results$IMMO == 
        "IM") == 6L && sum(patient_results$IMMO == "MO") == 2L)
    patient_results <- patient_results[order(patient_results$patient_id), , drop = FALSE]
    if (all(c("eligible_pair_n", "expected_pair_n") %in% names(patient_results))) 
        stopifnot(all(patient_results$eligible_pair_n == patient_results$expected_pair_n))
    v <- patient_results$patient_r
    stopifnot(all(is.finite(v)))
    im <- patient_results$IMMO == "IM"
    observed <- mean(v[im]) - mean(v[!im])
    allocations <- combn(8L, 2L)
    differences <- apply(allocations, 2, function(j) mean(v[-j]) - mean(v[j]))
    extreme <- abs(differences) >= abs(observed) - config$permutation_tolerance
    list(statistics = data.frame(IM_n = 6L, MO_n = 2L, IM_median = median(v[im]), MO_median = median(v[!im]), 
        mean_difference_IM_minus_MO = observed, exact_P = mean(extreme), allocations = ncol(allocations)))
}

# ---- Analysis ----
jobs <- read_tab(input("similarity_jobs.tsv")) # lineage, counts, metadata, hvg
config <- list(scale_factor = 10000, min_cells = 50, permutation_tolerance = 1e-12)
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  counts <- readRDS(job$counts) # corrected counts for this lineage, ALL RNA features
  meta <- read_tab(job$metadata) # cell_id, sample_id, patient_id, IMMO
  hvg <- read_tab(job$hvg)$gene # exact lineage-specific shared HVG2000, never whole-atlas HVGs
  profiles <- similarity_lesion_profiles(counts, meta, hvg, config)
  pairs <- similarity_within_patient_pairs(profiles)
  patients <- similarity_patient_summary(pairs, config)
  test <- similarity_exact_permutation(patients, config)
  write_tab(pairs, paste0(job$lineage, "_lesion_pairs.tsv"))
  write_tab(patients, paste0(job$lineage, "_patient_similarity.tsv"))
  write_tab(test$statistics, paste0(job$lineage, "_exact_test.tsv"))
}
# Source code supports these operations; final median/exact-test provenance is epithelial.
# Do not reinterpret every historical immune-lineage result as having used this later test.
