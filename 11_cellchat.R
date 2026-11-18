# Separate IM/MO CellChat networks on the exact 29-state population.
# Run from this folder: Rscript --vanilla 11_cellchat.R
source("00_paths.R")


# ---- Analysis ----
library(CellChat)
expr <- readRDS(input("communication_RNA_data.rds")) # log-normalized genes × cells
meta <- read_tab(input("communication_metadata.tsv")) # barcode, state, IMMO
states <- read_tab("resources/communication_states.tsv")
stopifnot(identical(colnames(expr), meta$barcode), setequal(meta$state, states$state),
          sum(states$lineage == "T") == 17, sum(states$lineage == "Mac") == 9, sum(states$lineage == "DC") == 3,
          all(table(factor(meta$state, levels = states$state), factor(meta$IMMO, levels = c("IM", "MO"))) >= 10))
for (group in c("IM", "MO")) {
  m <- meta[meta$IMMO == group, ]; rownames(m) <- m$barcode
  cc <- createCellChat(object = expr[, m$barcode, drop = FALSE], meta = m, group.by = "state")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc, variable.both = FALSE)
  cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE, population.size = FALSE,
                         nboot = 1000, seed.use = 1, Kh = 0.5, n = 1)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  write_tab(subsetCommunication(cc), paste0(group, "_CellChat_LR.tsv"))
  saveRDS(cc, output(paste0(group, "_CellChat.rds")))
}
# No PPI projection or downsampling. Within-network P values are not IM-versus-MO tests.
