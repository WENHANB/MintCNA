# Pointwise circular-shift test on observed local maxima (your MissInspect scorer)
# Returns cps + a small table. Uses BH-FDR (default) or Holm-FWER.
get_cps_pointwise <- function(sim, h = 15, penalty = "l2",
                              B = 400, control = c("FDR","FWER"),
                              q = 0.10, alpha = 0.05, seed = 1L) {
  control <- match.arg(control); set.seed(seed)
  n <- ncol(sim)
  
  # 1) Scores + observed candidate peaks (local maxima) -----------------------
  sc_obs <- multiMissInspect_scores(sim, h = h, penalty = penalty)
  peaks  <- sc_obs$local_maxima                    # data.frame: position, score
  if (!nrow(peaks)) return(list(cps = integer(0), table = data.frame()))
  pos <- peaks$position
  Hobs <- abs(peaks$score)
  m <- length(pos)
  
  # 2) Build null at the SAME positions via circular shifts -------------------
  circ <- function(n){ k <- sample.int(n,1)-1L; ((0:(n-1L)) + k) %% n + 1L }
  Hnull <- matrix(NA_real_, nrow = B, ncol = m)
  for (b in seq_len(B)) {
    idx <- circ(n)
    sc_b <- multiMissInspect_scores(sim[, idx, drop = FALSE], h = h, penalty = penalty)
    sb   <- sc_b$score                                # full score vector (length n)
    Hnull[b, ] <- abs(sb[pos])                        # evaluate at SAME positions
  }
  
  # 3) Pointwise permutation p-values + selection -----------------------------
  # add 1 in numerator/denominator for unbiased permutation p-values
  p <- (1 + colSums(Hnull >= matrix(Hobs, B, m, byrow = TRUE))) / (B + 1)
  
  if (control == "FDR") {
    padj <- p.adjust(p, "BH"); sel <- which(padj <= q)
  } else { # FWER
    padj <- p.adjust(p, "holm"); sel <- which(padj <= alpha)
  }
  
  cps <- sort(pos[sel])
  tab <- data.frame(t = pos[sel], height = Hobs[sel],
                    p = p[sel], p_adj = padj[sel],
                    method = if (control=="FDR") "pointwise_BH" else "pointwise_Holm")
  list(cps = cps, table = tab)
}
