# =====================================================================
# R6  Evidence Hierarchy and Inferential Decision Rules     (Section 3.4.6)
# ---------------------------------------------------------------------
# Numbered equations: none. Section 3.4.6 states an evidence hierarchy and five
# decision rules; this script encodes them as explicit, auditable code that reads
# the outputs of R3, R4 and R5 and attaches the wording each result is allowed to carry.
#
#   Tier 1  Principal associational evidence : baseline Dg coefficient (3.12), post-2020
#           interaction (3.14); conventional city-clustered p AND wild-bootstrap p.
#           Where they disagree, the bootstrap decides statistical certainty.
#   Tier 2  Descriptive dynamic evidence      : Moran/LISA (R4), Markov (R5): no causal language.
#   Tier 3  Measurement/specification sensitivity : alt. HQD, LOCO, LODO, trends, FD,
#           strength-balance: report DIRECTIONAL consistency, never "robustness".
#   Tier 4  Secondary conditional evidence    : industrial interactions, TP-only, SLX,
#           HSR thresholds (raw, Holm, BH).
#   Rules 1-5 are implemented in classify_inference(), directional_consistency(),
#   classify_secondary(), and the fixed wording column.
# =====================================================================
setwd("C:/R/spatial") 
source("7. R_code for reproducibility.R")
need_output(tag("results_R3.rds"), "R3_panel_regression_sensitivity.R")
res <- readRDS(outp(tag("results_R3.rds")))
ALPHA10 <- 0.10

# ---- Rules 1, 2, 5: classification of a principal coefficient --------------------------------
classify_inference <- function(est, p_conv, p_boot, alpha = ALPHA10) {
  conv <- p_conv < alpha; boot <- p_boot < alpha
  if (conv && boot)  "Statistically supported under both conventional (city-clustered) and wild-bootstrap inference"
  else if (conv)     "Weakly/suggestively supported under conventional inference and inferentially fragile overall (not confirmed by wild bootstrap); not robustly significant"
  else if (boot)     "Supported by wild bootstrap only; treat conventional-inference discrepancy with caution"
  else               "Insufficient evidence to distinguish the coefficient from zero (not proof of absence or equality)"
}
# ---- Rule 2: directional consistency across sensitivity models ---------------------------------
directional_consistency <- function(benchmark_est, est_vec) {
  est_vec <- est_vec[!is.na(est_vec)]
  data.frame(n_models = length(est_vec), n_same_sign = sum(sign(est_vec) == sign(benchmark_est)),
             share_same_sign = mean(sign(est_vec) == sign(benchmark_est)),
             wording = "Directional consistency only; not statistical robustness")
}
# ---- Secondary evidence (Rule 4 + multiple-testing audit) -------------------------------------------
classify_secondary <- function(p_raw, p_holm, p_bh, alpha = ALPHA10) {
  if (is.na(p_raw)) return("Not estimable in the available sample")
  if (!is.na(p_holm) && p_holm < alpha) "Survives Holm (FWER) adjustment; still secondary associational evidence"
  else if (!is.na(p_bh) && p_bh < alpha) "Survives Benjamini-Hochberg (FDR) but not Holm; secondary associational evidence"
  else if (p_raw < alpha) "Raw significance only; does not survive multiple-testing adjustment"
  else "No stable evidence under the specified definition (not proof that cities do not interact)"
}

# ---- Tier 1 -------------------------------------------------------------------------------------------
boot <- read.csv(outp(tag("results_main_inference.csv")))
tier1 <- data.frame(tier = "1 Principal associational", item = boot$model, estimate = boot$est,
                    p_conventional = boot$p_conv, p_bootstrap = boot$p_boot,
                    ci95 = sprintf("[%.4f, %.4f]", boot$ci95_lo, boot$ci95_hi),
                    verdict = mapply(classify_inference, boot$est, boot$p_conv, boot$p_boot))

# ---- Tier 3 (directional consistency vs benchmark sign) -----------------------------------------------------
bench <- res$base$est[res$base$term == "Dg"]
sens_blocks <- list(`alternative HQD (baseline rows)` = res$alt[grepl("^baseline", res$alt$model) & res$alt$term == "Dg", "est"],
                    `leave-one-city-out` = res$loco$est[grepl("^HQD_ent", res$loco$model)],
                    `leave-one-dimension-out` = res$lodo$est,
                    `city trends (3.15)` = res$trend$est, `first differences (3.16)` = res$fd$est)
tier3 <- do.call(rbind, Map(function(v, nm) cbind(tier = "3 Sensitivity", item = nm, directional_consistency(bench, v)),
                            sens_blocks, names(sens_blocks)))

# ---- Tier 4 -----------------------------------------------------------------------------------------------------
# script 3 version (e) labels the models "HQD_ent 3.25 ..." / "HQD_eq 3.25 ..."
cond <- res$cond[grepl("3.25", res$cond$model, fixed = TRUE) & grepl("^INT_", res$cond$term), ]
tp   <- res$cond[grepl("3.30", res$cond$model, fixed = TRUE) & grepl("^INT_", res$cond$term), ]
tier4 <- data.frame(tier = "4 Secondary conditional", item = c(cond$model, tp$model), estimate = c(cond$est, tp$est),
                    p_conventional = c(cond$p, tp$p),
                    verdict = c(ifelse(cond$p < ALPHA10, "Conventional significance only (no bootstrap run for secondary models); conditional heterogeneity, not causal moderation",
                                       "No sufficient evidence of interaction"),
                                ifelse(tp$p < ALPHA10, "Interaction persists without industry subsystem (reduces construct-overlap concern)",
                                       "Interaction does not persist without industry subsystem")))
thr_file <- c(outp("results_hsr_thresholds.csv"), outp("DEMO_results_hsr_thresholds.csv"))
thr_file <- thr_file[file.exists(thr_file)][1]
if (!is.na(thr_file)) {
  thr <- read.csv(thr_file)
  tier4 <- rbind(tier4[, c("tier", "item", "estimate", "p_conventional", "verdict")],
                 data.frame(tier = "4 Secondary conditional", item = paste("HSR", thr$outcome, thr$cutoff, thr$side),
                            estimate = thr$est, p_conventional = thr$p,
                            verdict = mapply(classify_secondary, thr$p, thr$p_holm, thr$p_bh)))
}

# ---- Rules 3 and 4: descriptive evidence carries a fixed wording ----------------------------------------------------
tier2 <- data.frame(tier = "2 Descriptive dynamic", item = c("Global Moran's I / LISA (R4)", "TPI, HQD, Joint Markov (R5)"),
                    verdict = "Descriptive only: spatial/transition patterns; high diagonal = persistence, not lock-in; no causal language")

pad <- function(x, cols) { for (c in setdiff(cols, names(x))) x[[c]] <- NA; x[, cols] }
cols <- c("tier", "item", "estimate", "p_conventional", "p_bootstrap", "ci95", "n_models", "n_same_sign", "share_same_sign", "verdict", "wording")
tier3$verdict <- tier3$wording
ev <- rbind(pad(tier1, cols), pad(tier2, cols), pad(tier3, cols), pad(tier4, cols))
write.csv(ev, outp(tag("evidence_hierarchy.csv")), row.names = FALSE)
cat("\nEvidence hierarchy (decision rules applied):\n"); print(ev[, c("tier", "item", "estimate", "p_conventional", "p_bootstrap", "verdict")], right = FALSE, row.names = FALSE)
