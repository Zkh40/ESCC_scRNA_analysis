# ESCC scRNA-seq core analysis

Numbered analysis scripts for the multifocal ESCC study (8 patients, 17 lesions; IM 12, MO 5). Start here. Parameters are visible in the scripts; `KEY_PARAMETERS.tsv` is the short reading/checking index.

## Workflow

| Script | Input | Main output |
|---|---|---|
| 01_qc_doublets.R | Original 10x UMI matrices; matrix_manifest.tsv | qc_singlets.rds |
| 02_integrate_annotate.R | 01 + reviewed broad labels | atlas_annotated.rds; cluster table/UMAP |
| 03_immune_subclusters.R | 02 + reviewed immune labels | TNK_clean.rds; Myeloid_clean.rds |
| 04_signature_scores.R | Exact TNK/Mac/DC scoring parents; supplied gene panels | Scores; complete 60-test Ttex/Tpex/TRM comparison |
| 05_de_kegg.R | Selected T/DC/Mac populations; original feature universe; frozen KEGG tables | DE and ORA/GSEA tables |
| 06_cluster_expression_gsea.R | Frozen TNK expression population and requested cluster/signature pairs | Pooled-condition expression enrichment |
| 07_epithelial_copykat.R | Original UMI; same-sample epithelial queries/TNK references | CopyKAT native files and predictions |
| 08_epithelial_state.R | Reviewed frozen malignant epithelial object; Hallmark50 | t-SNE, scores, CytoTRACE, ITH, GSVA |
| 09_cd8_trajectory.R | CD8 count/labels OR the historical tree | New estimate OR historical rerooted pseudotime |
| 10_abundance.R | Two explicit lesion/state count tables | Fractions and Spearman correlations |
| 11_cellchat.R | Frozen 29-state RNA/data and metadata | Separate IM/MO networks and LR tables |
| 12_transcript_similarity.R | Each lineage's corrected counts, metadata and frozen shared HVG2000 | Lesion pairs, patient medians, exact comparison |

01→02→03 is the upstream sequence. 04–12 are analysis branches with their own input populations, not a single automatic chain. The historical malignant-like subset includes reviewed identity evidence; CopyKAT predictions alone do not recreate that selection. This release starts with expression matrices, not FASTQ processing.

## Running

Use the existing study environment in `software_versions.tsv` (main: R 4.4.2 / Seurat 5.3.1). Run each branch in a new R process. Monocle2 requires its isolated dplyr 0.8.5 profile; CellChat and CytoTRACE require their recorded profiles. No installation or automatic environment switching occurs here.

1. Set `data_dir` and a **new, unused** `results_dir` in `00_paths.R` (or ESCC_DATA / ESCC_RESULTS).
2. Supply the input files below. Run from this package's root, for example `Rscript --vanilla 01_qc_doublets.R`.
3. Review the exported cluster/marker evidence before supplying `reviewed_broad_labels.tsv` or `reviewed_immune_labels.tsv`. Supplied resource labels are historical references, not automatic labels for a newly clustered object.
4. Run the required downstream branch with its original frozen population. Scripts save ordinary RDS/TSV/PDF outputs. Never point results_dir at the frozen project outputs or run twice into an existing result folder.

For trajectory, default mode is historical rerooting. Set ESCC_TRAJECTORY_MODE=estimate to estimate a new tree instead. Never apply historical State3/Y_78 to a newly estimated tree.

## External input schemas

All TSVs have headers. Object paths, matrix paths and feature-list paths inside TSV files must be absolute. Expression matrices are genes × cells with unique names. Seurat objects use RNA `counts` and `data`; `data` means log-normalized expression. Metadata is aligned by cell ID, not row position guessed from another file. Gene lists have a `gene` column, except original-feature lists with `SYMBOL`.

| File in data_dir | Required contents |
|---|---|
| matrix_manifest.tsv | sample_id, matrix_dir; all 17 sample IDs in resources/samples.tsv, absolute 10x directories |
| reviewed_broad_labels.tsv | cluster, label; reviewed global-cluster mapping |
| reviewed_immune_labels.tsv | lineage (TNK/Myeloid), cluster (C-number), label, retained (TRUE/FALSE); include rejected clusters |
| TNK_score_parent.rds / Mac_score_parent.rds / DC_score_parent.rds | Original scoring parents, metadata cluster and IMMO; respectively clean TNK, Mac14365, and DC C13 753 cells |
| de_jobs.tsv | id, profile (T/DC/Mac), object, feature_symbols; object is exact selected Seurat subset; feature_symbols is full original parent RNA feature TSV |
| TERM2GENE_hsa.tsv / TERM2NAME_hsa.tsv | Frozen 2026-09-11 KEGG snapshot, columns term/entrez and term/name; supply authorized copies |
| <DE job id>_symbol_entrez.tsv | Mac jobs only; original SYMBOL/ENTREZID mapping, including ambiguous alternatives |
| TNK_expression_parent.rds / expression_gsea_jobs.tsv | Frozen TNK Seurat object (cluster, IMMO); jobs: id, cluster, signature |
| copykat_cells.tsv | sample_id, barcode_raw, role (query/reference); references are T/NK from the same sample |
| epithelial_malignant.rds / hallmark50.rds | Frozen malignant-like Seurat object with copykat_malignancy_call=MALIGNANT_LIKE, sample_id, IMMO and PCA; named list of 50 frozen Hallmark sets |
| cd8_counts.rds / cd8_labels.tsv | Estimation only: RNA counts; barcode, cluster in exact column order, seven declared CD8 clusters |
| cd8_historical_tree.rds | Rerooting only: original Monocle2 tree with State3/Y_78 |
| global_abundance_counts.tsv / targeted_abundance_counts.tsv | sample_id, state, count, denominator; complete 17×38 and 17×2 grids, including zeros; denominators as documented in script10 |
| communication_RNA_data.rds / communication_metadata.tsv | Frozen genes×cells logRNA matrix; barcode, state, IMMO; exact 29 states in resources/communication_states.tsv |
| similarity_jobs.tsv | lineage, counts, metadata, hvg; counts RDS is corrected counts for all features; metadata TSV has cell_id/sample_id/patient_id/IMMO; hvg TSV has 2000 frozen genes |

Study data/Seurat objects, the external KEGG snapshot, full Hallmark50 and the lineage-specific shared HVG2000 lists are not redistributed in this code archive. Input IDs and population definitions must match the original analysis. No objects are inferred from filenames or downloaded automatically.
