#!/usr/bin/env Rscript
# CDKN2A/B and CCND2 at single-cell resolution, against the matched bulk WGS.
#
# Johnson et al. (Nature 2026), the CARE cohort's own paper, report that acquired
# CDKN2A deletion and CCND2 amplification at recurrence track with an increase in
# proliferating stem-like cells, and that recurrence is BIMODAL -- some tumours
# stay genetically stable, others acquire new alterations. Both loci are named
# there, so finding them independently is a check on the method rather than a
# fishing expedition.
suppressPackageStartupMessages(library(data.table))
setwd("<MINTCNA_ROOT>/ATAC_RNA")
TAG <- "h15a010e020"
LOCI <- data.table(
  gene = c("CDKN2A", "CDKN2B", "CCND2"),
  chr  = c("chr9", "chr9", "chr12"),
  from = c(21967752, 22002903, 4273761),      # hg38
  to   = c(21995301, 22009313, 4305353),
  dir  = c("loss", "loss", "gain"))
S <- c("P59T2","P59T3","P100T1","P100T3","P103T1","P103T2","P104T1","P104T2","P105T1","P105T3",
       "P106T1","P106T2","P108T1","P108T2","P109T1","P109T2","P110T1","P110T2","P111T2","P111T3")
man <- fread("Data/meta/sample_manifest.csv"); setkey(man, care_id)
ORD <- c("Initial"=1, "Recurrence"=2, "2nd Recurrence"=3)

res <- rbindlist(lapply(S, function(s) {
  e <- new.env(); load(sprintf("Result/MintCNA/%s_multi/%s_results_%s.RData", s, s, TAG), e)
  Sm <- base::get("cna_state_matrix", e); rf <- as.data.table(base::get("ref", e))
  keep <- fread(sprintf("Data/cell_lists/%s_cells.tsv", s))[cell_type=="Malignant", barcode]
  i <- colnames(Sm) %in% keep
  tr <- fread(sprintf("Result/truth/%s_truth_100kb.tsv.gz", s))
  rbindlist(lapply(seq_len(nrow(LOCI)), function(k) {
    L <- LOCI[k]
    # every MintCNA bin overlapping the gene body
    b <- which(rf$chr == L$chr & rf$end >= L$from & rf$start <= L$to)
    # the truth's own 100 kb bins over the same span
    t <- tr[chr == L$chr & end >= L$from & start <= L$to]
    if (!length(b)) return(NULL)
    st <- Sm[b, i, drop=FALSE]
    data.table(care_id = s, patient = man[s, patient], tp = man[s, timepoint],
               ord = ORD[man[s, timepoint]], subtype = substr(man[s, subtype], 1, 5),
               gene = L$gene, n_cell = sum(i), n_bin = length(b),
               # fraction of malignant nuclei calling the event in ANY overlapping bin
               frac = if (L$dir == "loss") mean(colSums(st == 1) > 0) else mean(colSums(st == 3) > 0),
               wgs_log2 = mean(t$log2, na.rm = TRUE),
               wgs_call = paste(unique(t$cnv_call[t$cnv_call %in% c("-","0","+")]), collapse="/"))
  }))
}))
fwrite(res, "Result/qc/locus_cdkn2a_ccnd2.tsv", sep="\t")

for (g in c("CDKN2A","CCND2")) {
  cat(sprintf("\n=== %s (%s) ===\n\n", g, LOCI[gene==g, dir]))
  d <- res[gene == g]; setorder(d, patient, ord)
  print(d[, .(patient, care_id, tp, subtype, n_bin,
              pct_cells = round(100*frac, 1),
              wgs_log2 = round(wgs_log2, 3), wgs_call)], nrows=25)
  p <- d[, .(early = care_id[1], late = care_id[.N], subtype = subtype[1],
             f_e = frac[1], f_l = frac[.N], w_e = wgs_log2[1], w_l = wgs_log2[.N]), by = patient]
  cat("\n  paired, earlier -> later:\n")
  print(p[, .(patient, subtype, pct_e = round(100*f_e,1), pct_l = round(100*f_l,1),
              d_pct = round(100*(f_l-f_e),1), wgs_e = round(w_e,3), wgs_l = round(w_l,3))], nrows=12)
}
