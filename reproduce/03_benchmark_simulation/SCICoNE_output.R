library(data.table)

args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else "c3_N100_100kb_l1e4"
print(BASE)

CNA_DIR <- file.path("<MINTCNA_ROOT>/Simulation/cnasim", BASE)

cna_scicone <- t(fread(file.path(CNA_DIR,"Calling/SCICoNE/run1_inferred_cnvs.csv"), header = F))

cna_scicone <- cna_scicone + 1

colnames(cna_scicone) <- paste0("cell", 1:ncol(cna_scicone))

save(cna_scicone, file = file.path(CNA_DIR,"Calling/SCICoNE.RData"))

table(cna_scicone)
