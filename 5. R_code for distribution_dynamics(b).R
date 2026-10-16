# =====================================================================
# R5  Distribution-dynamics Analysis (Markov chains)       (Section 3.4.5)
#     VERSION (b): prints Tables 4.3, 4.4, 5.2, 5.3, 5.4, 5.5, A2A, A2B, A2C
#     and compares every number with the thesis (column "Same_as_thesis").
# ---------------------------------------------------------------------
# Equations implemented
#   (3.48) pooled quartile states (4 states) from the full 2011-2023 distribution;
#          pooled terciles (3 states) for the Joint Markov analysis
#   (3.49) first-order transition probability p_rs = n_rs / sum_k n_rk
#   (3.50) joint state  J_it = 3 (d_it - 1) + h_it,  J in {1,...,9}
# What changed compared with version (a)
#   1. Cut-points use the MIDPOINT rule (q_mid below). This reproduces the
#      thesis values exactly (Table 4.3, A2A). R's quantile(type = 7) gives
#      slightly different cut-points.
#   2. The cut-points are kept and printed (Tables 4.3, 5.2, A2A).
#   3. summ_rates() adds unweighted persistence, adjacent-up and adjacent-down
#      shares (Tables 4.4 and 5.3).
#   4. The HQD and Joint transition matrices are printed with counts
#      (Figures 5.1, 5.2; Tables A2B, A2C).
#   5. conv() adds p99 and N9 (Table 5.5); the joint-state map is printed (Table 5.4).
# TPI Markov uses the real data (Dg). HQD and Joint Markov need the HQD index:
# real if hqd_controls.xlsx is supplied, otherwise DEMO_ output.
# =====================================================================
setwd("C:/R/spatial")
source("7. R_code for reproducibility.R")
options(scipen = 999)
need_output("tpi_index.rds", "R2_coupling_coordination.R")
need_output(tag("hqd_indices.rds"), "R1_standardisation_entropy.R")
dat <- merge(readRDS(outp("tpi_index.rds"))[, c("city_code", "year", "Dg")],
             readRDS(outp(tag("hqd_indices.rds")))[, c("city_code", "year", "HQD_ent")], by = c("city_code", "year"))
dat <- dat[order(dat$city_code, dat$year), ]

# ---- (3.48) pooled cut-points and states --------------------------------------
# Midpoint quantile: sort the 208 values, find position (n - 1) * p, and take
# the average of the two values on either side of that position.
q_mid <- function(z, p) {
  s   <- sort(z)
  pos <- (length(s) - 1) * p
  (s[floor(pos) + 1] + s[ceiling(pos) + 1]) / 2
}
cut_points <- function(z, K) {                 # K - 1 cut-points: quartiles (K = 4), terciles (K = 3)
  q <- numeric(K - 1)
  for (k in 1:(K - 1)) q[k] <- q_mid(z, k / K)
  q
}
state_from_cuts <- function(z, q) {            # state 1 if z <= q1, ..., state K if z > q_{K-1}
  as.integer(1 + rowSums(outer(z, q, ">")))
}

CUT <- list(tpi4 = cut_points(dat$Dg, 4),           # Table 4.3
            hqd4 = cut_points(dat$HQD_ent, 4),      # Table 5.2
            tpi3 = cut_points(dat$Dg, 3),           # Table A2A
            hqd3 = cut_points(dat$HQD_ent, 3))      # Table A2A

dat$S_tpi4 <- state_from_cuts(dat$Dg,      CUT$tpi4)
dat$S_hqd4 <- state_from_cuts(dat$HQD_ent, CUT$hqd4)
dat$d3     <- state_from_cuts(dat$Dg,      CUT$tpi3)
dat$h3     <- state_from_cuts(dat$HQD_ent, CUT$hqd3)
dat$J      <- 3 * (dat$d3 - 1) + dat$h3                              # (3.50)
stopifnot(all(dat$J %in% 1:9))

# ---- consecutive-year pairs (each city contributes 12 transitions) ---------------
pairs <- merge(dat, transform(dat, year = year - 1), by = c("city_code", "year"), suffixes = c("", "_next"))
stopifnot(nrow(pairs) == 16 * 12)
pairs$period <- ifelse(pairs$year <= 2019, "origin 2011-2019", "origin 2020-2022")   # 2019->2020 belongs to the earlier group
PERIODS <- c("full period", "origin 2011-2019", "origin 2020-2022")

# ---- (3.49) transition matrices and summary rates -------------------------------
trans_counts <- function(from, to, K) table(factor(from, 1:K), factor(to, 1:K))
trans_prob   <- function(N) { rs <- rowSums(N); N / ifelse(rs == 0, NA, rs) }      # (3.49)

summ_rates <- function(N) {
  N  <- matrix(as.vector(N), nrow = nrow(N))     # plain count matrix
  r  <- row(N); c <- col(N); tot <- sum(N); rs <- rowSums(N)
  diag_p <- diag(N)[rs > 0] / rs[rs > 0]          # staying probability of each origin state that occurs
  data.frame(transitions            = tot,
             persistence_unweighted = mean(diag_p),                 # simple average of the diagonal
             persistence            = sum(N[r == c]) / tot,         # weighted by origin counts
             upward                 = sum(N[c > r]) / tot,
             downward               = sum(N[c < r]) / tot,
             adjacent_up            = sum(N[c == r + 1]) / tot,     # up by exactly one state
             adjacent_down          = sum(N[c == r - 1]) / tot,     # down by exactly one state
             skip_move              = sum(N[abs(c - r) > 1]) / tot)
}
run_markov <- function(from, to, K, label, per) {
  out <- list(); S <- list()
  for (p in PERIODS) {
    i <- if (p == "full period") rep(TRUE, length(from)) else per == p
    N <- trans_counts(from[i], to[i], K); P <- trans_prob(N)
    out[[p]] <- list(N = N, P = P); S[[p]] <- cbind(variable = label, period = p, summ_rates(N))
  }
  list(mats = out, summary = do.call(rbind, S))
}
mk_long <- function(res, label) do.call(rbind, lapply(names(res$mats), function(p) {
  N <- res$mats[[p]]$N; P <- res$mats[[p]]$P
  data.frame(variable = label, period = p, from = rep(rownames(N), ncol(N)), to = rep(colnames(N), each = nrow(N)),
             n = as.vector(N), n_origin = rep(rowSums(N), ncol(N)), prob = as.vector(P)) }))

mk_tpi   <- run_markov(pairs$S_tpi4, pairs$S_tpi4_next, 4, "TPI (Dg), quartiles",       pairs$period)
mk_hqd   <- run_markov(pairs$S_hqd4, pairs$S_hqd4_next, 4, "HQD (entropy), quartiles",  pairs$period)
mk_joint <- run_markov(pairs$J,      pairs$J_next,      9, "Joint TPI x HQD, terciles", pairs$period)

# ---- (3.50) Joint Markov conversion metrics (Table 5.5) ---------------------------------
dd <- pairs$d3_next - pairs$d3; dh <- pairs$h3_next - pairs$h3
conv <- function(i) {
  N9 <- sum(pairs$J[i] == 9)                                       # transitions starting in joint state 9
  data.frame(transitions              = sum(i),
             same_joint_state         = mean(pairs$J[i] == pairs$J_next[i]),
             adjacent_synchronised_up = mean(dd[i] == 1 & dh[i] == 1),   # (d,h) -> (d+1,h+1)
             tpi_only_adjacent_up     = mean(dd[i] == 1 & dh[i] == 0),   # (d,h) -> (d+1,h)
             hqd_only_adjacent_up     = mean(dh[i] == 1 & dd[i] == 0),   # (d,h) -> (d,h+1)
             p99 = if (N9 > 0) sum(pairs$J[i] == 9 & pairs$J_next[i] == 9) / N9 else NA,
             N9  = N9,
             any_up_tpi = mean(dd[i] > 0), any_up_hqd = mean(dh[i] > 0))
}
conv_tab <- NULL
for (p in PERIODS) {
  i <- if (p == "full period") rep(TRUE, nrow(pairs)) else pairs$period == p
  conv_tab <- rbind(conv_tab, cbind(period = p, conv(i)))
}

# =====================================================================
# THESIS VALUES (used only for the comparison columns)
# =====================================================================
T_CUTS <- data.frame(
  table  = c(rep("Table 4.3  TPI quartiles", 3), rep("Table 5.2  HQD quartiles", 3),
             rep("Table A2A  TPI terciles", 2),  rep("Table A2A  HQD terciles", 2)),
  cut    = c("q25", "q50", "q75", "q25", "q50", "q75", "q33", "q67", "q33", "q67"),
  R      = c(CUT$tpi4, CUT$hqd4, CUT$tpi3, CUT$hqd3),
  thesis = c(0.483491, 0.535048, 0.604794, 0.094188, 0.183604, 0.221444,
             0.504852, 0.574045, 0.159145, 0.205923))

RATE_COLS <- c("transitions", "persistence_unweighted", "persistence", "upward", "downward",
               "adjacent_up", "adjacent_down")
T44 <- data.frame(period = PERIODS, transitions = c(192, 144, 48),          # Table 4.4 (TPI)
                  persistence_unweighted = c(0.738, 0.757, 0.506), persistence = c(0.734, 0.743, 0.708),
                  upward = c(0.188, 0.194, 0.167), downward = c(0.078, 0.062, 0.125),
                  adjacent_up = c(0.172, 0.174, 0.167), adjacent_down = c(0.078, 0.062, 0.125))
T53 <- data.frame(period = PERIODS, transitions = c(192, 144, 48),          # Table 5.3 (HQD)
                  persistence_unweighted = c(0.789, 0.764, 0.840), persistence = c(0.786, 0.757, 0.875),
                  upward = c(0.177, 0.222, 0.042), downward = c(0.036, 0.021, 0.083),
                  adjacent_up = c(0.167, 0.208, 0.042), adjacent_down = c(0.031, 0.014, 0.083))
CONV_COLS <- c("transitions", "adjacent_synchronised_up", "tpi_only_adjacent_up",
               "hqd_only_adjacent_up", "p99", "N9")
T55 <- data.frame(period = PERIODS, transitions = c(192, 144, 48),          # Table 5.5
                  adjacent_synchronised_up = c(0.010, 0.014, 0),
                  tpi_only_adjacent_up = c(0.099, 0.104, 0.083),
                  hqd_only_adjacent_up = c(0.109, 0.146, 0),
                  p99 = c(0.906, 0.947, 0.846), N9 = c(32, 19, 13))

# "yes" when every listed column agrees with the thesis to its printed precision
same_rows <- function(R, thesis, cols, tol = 0.0006) {
  ok <- rep(TRUE, nrow(R))
  for (cl in cols) ok <- ok & !is.na(R[[cl]]) & abs(R[[cl]] - thesis[[cl]]) <= tol
  ifelse(ok, "yes", "NO")
}
T_CUTS$Same_as_thesis <- ifelse(abs(T_CUTS$R - T_CUTS$thesis) <= 0.0000015, "yes", "NO")

# Transition matrix shown as "count [probability]" plus origin totals (Tables A2B, A2C)
count_prob <- function(N) {
  K  <- nrow(N)
  N  <- matrix(as.numeric(as.vector(N)), nrow = K)
  rs <- rowSums(N)
  P  <- N / ifelse(rs == 0, NA, rs)
  cell <- ifelse(is.na(P), paste0(N, " [-]"), sprintf("%.0f [%.3f]", N, P))
  out <- data.frame(matrix(cell, nrow = K), origin_N = rs)
  names(out) <- c(paste0("to ", 1:K), "origin_N")
  rownames(out) <- paste0("from ", 1:K)
  out
}

# =====================================================================
# PRINT
# =====================================================================
if (DEMO_EXTRA) cat("\nNOTE: hqd_controls.xlsx not found, so HQD is SIMULATED.\n",
                    "     Only the TPI rows (Tables 4.3, 4.4, A2A TPI, A2B TPI) can be compared with the thesis.\n")

cat("\n========== Tables 4.3, 5.2, A2A: pooled cut-points (3.48) ==========\n")
print(data.frame(T_CUTS[, 1:2], R = round(T_CUTS$R, 6), thesis = T_CUTS$thesis,
                 Same_as_thesis = T_CUTS$Same_as_thesis), row.names = FALSE)

cat("\n========== Table A2B: TPI transitions, count [probability] (3.49) ==========\n")
for (p in PERIODS) { cat("\n--", p, "--\n"); print(count_prob(mk_tpi$mats[[p]]$N)) }

cat("\n========== Table 4.4: TPI Markov summary ==========\n")
s_tpi <- mk_tpi$summary
s_tpi$Same_as_thesis <- same_rows(s_tpi, T44, RATE_COLS)
print(s_tpi[, c("period", RATE_COLS, "Same_as_thesis")], digits = 3, row.names = FALSE)

cat("\n========== Figure 5.1 / Table A2B: HQD transitions, count [probability] (3.49) ==========\n")
for (p in PERIODS) { cat("\n--", p, "--\n"); print(count_prob(mk_hqd$mats[[p]]$N)) }

cat("\n========== Table 5.3: HQD Markov summary ==========\n")
s_hqd <- mk_hqd$summary
s_hqd$Same_as_thesis <- same_rows(s_hqd, T53, RATE_COLS)
print(s_hqd[, c("period", RATE_COLS, "Same_as_thesis")], digits = 3, row.names = FALSE)
cat("Thesis Table 5.3:\n"); print(T53, row.names = FALSE)

cat("\n========== Table 5.4: joint state map J = 3(d - 1) + h (3.50) ==========\n")
map54 <- data.frame(joint_state = 1:9, D_state = rep(1:3, each = 3), HQD_state = rep(1:3, times = 3))
stopifnot(all(map54$joint_state == 3 * (map54$D_state - 1) + map54$HQD_state))
print(map54, row.names = FALSE)
cat("City-years in each joint state:\n"); print(table(factor(dat$J, 1:9)))

cat("\n========== Figure 5.2 / Table A2C: Joint transitions, full period, count [probability] ==========\n")
print(count_prob(mk_joint$mats[["full period"]]$N))

cat("\n========== Table 5.5: Joint Markov metrics (3.50) ==========\n")
conv_tab$Same_as_thesis <- same_rows(conv_tab, T55, CONV_COLS)
print(conv_tab[, c("period", CONV_COLS, "Same_as_thesis")], digits = 3, row.names = FALSE)
cat("Thesis Table 5.5:\n"); print(T55, row.names = FALSE)

# =====================================================================
# SAVE
# =====================================================================
write.csv(T_CUTS, outp(tag("table_4_3_5_2_A2A_cut_points.csv")), row.names = FALSE)
write.csv(mk_long(mk_tpi, "TPI"), outp("markov_tpi_transitions.csv"), row.names = FALSE)           # TPI = real data
write.csv(mk_long(mk_hqd, "HQD"), outp(tag("markov_hqd_transitions.csv")), row.names = FALSE)
write.csv(mk_long(mk_joint, "Joint"), outp(tag("markov_joint_transitions.csv")), row.names = FALSE)
write.csv(s_tpi, outp("table_4_4_tpi_markov_summary.csv"), row.names = FALSE)
write.csv(s_hqd, outp(tag("table_5_3_hqd_markov_summary.csv")), row.names = FALSE)
write.csv(map54, outp("table_5_4_joint_state_map.csv"), row.names = FALSE)
write.csv(conv_tab, outp(tag("table_5_5_joint_markov_metrics.csv")), row.names = FALSE)
write.csv(rbind(mk_tpi$summary, mk_hqd$summary, mk_joint$summary), outp(tag("markov_summary.csv")), row.names = FALSE)
write.csv(dat, outp(tag("markov_states.csv")), row.names = FALSE)
cat("\nTables saved in:", normalizePath(CFG$out_dir), "\n")
