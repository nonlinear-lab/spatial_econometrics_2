# =====================================================================
# R3  Panel Regression, Temporal Heterogeneity, Sensitivity Analysis
#     and spatial weight matrices: Equations (3.12) - (3.37)
# ---------------------------------------------------------------------
# COMPLETE SCRIPT: run it as it is with source(). Nothing to paste.
#   Part A  setup and data (your original lines)
#   Part B  regressions written as plain lm() calls; matrices in base R
#   Part C  print every result to the console and save it as CSV
# Needs once:  install.packages("sandwich")
# =====================================================================

# =====================================================================
# PART A  SETUP AND DATA
# =====================================================================
setwd("C:/R/spatial") 
source("7. R_code for reproducibility.R")
need_output("tpi_index.rds", "R2_coupling_coordination.R")
need_output(tag("hqd_indices.rds"), "R1_standardisation_entropy.R")

tp <- readRDS(outp("tpi_index.rds")); hq <- readRDS(outp(tag("hqd_indices.rds")))
hq <- hq[, setdiff(names(hq), c("I1", "I2"))]     # I1, I2 are taken from the TPI panel below (avoids I1.x / I1.y)
raw <- readRDS(outp("tpi_panel.rds"))[, c("city_code", "year", "I1", "I2")]
dat <- merge(merge(tp, hq, by = c("city_code", "year")), raw, by = c("city_code", "year")); dat <- dat[order(dat$city_code, dat$year), ]
stopifnot(nrow(dat) == 208)

# (Estimation engine -- design(), fit_cl(), twfe(), coef_tab(), wcr_p(), wcr_ci(), boot_row() -- lives in
#  R7_reproducibility.R so that R3 and R4 share exactly the same computational implementation.)

# ======================= variables ===============================
dat$Dg_post <- dat$Dg * as.numeric(dat$year >= 2020)             # D_g x Post2020
save_csv <- function(x, name) { write.csv(x, outp(tag(name)), row.names = FALSE); invisible(x) }
res <- list()

# =====================================================================
# PART B  REGRESSIONS (3.12) - (3.30), each one a plain lm() call
# ---------------------------------------------------------------------
# Two-way fixed effects are written directly in the lm() formula:
#   mu_i     -> factor(city_code)   (city dummies)
#   lambda_t -> factor(year)        (year dummies)
# The slope estimates are identical to the fixed-effects (within) estimator.
# lm() does not cluster standard errors, so report() below takes each
# fitted lm and produces city-clustered SEs (G = 16), p-values from t(G-1),
# and the within R-squared.
# =====================================================================
library(sandwich)

CTRL <- c("Con1", "Con2", "Con3", "Con4", "Con5")   # names used in the result tables

# ---- clustered covariance matrix of a fitted lm ---------------------
# CR1: G/(G-1) * (N-1)/(N-K). With adjust = "fe_nested" the city dummies are
# not counted in K (they are nested in the city clusters), which is the
# Stata xtreg/reghdfe rule. Use adjust = "lm" to count them.
vcov_city <- function(m, adjust = "fe_nested") {
  cl <- expand.model.frame(m, ~ city_code)$city_code
  N  <- nobs(m); G <- length(unique(cl)); K <- m$rank
  V  <- vcovCL(m, cluster = cl, type = "HC1")
  if (adjust == "fe_nested" && "factor(city_code)" %in% attr(terms(m), "term.labels"))
    V <- V * (N - K) / (N - (K - (G - 1)))
  list(V = V, N = N, G = G)
}

# ---- result table for a fitted lm ------------------------------------
report <- function(m, label, terms) {
  vc <- vcov_city(m)
  b  <- coef(m)[terms]
  se <- sqrt(diag(vc$V))[terms]
  t  <- b / se
  p  <- 2 * pt(-abs(t), df = vc$G - 1)
  # within R-squared, two definitions:
  #  r2_within : removes city AND year effects (and city trends) -> like Table 5.13 (about 0.2)
  #  r2_city   : removes city effects only; year dummies count as explained variation
  #              -> the "Within R2" printed in Tables 5.6, 5.7, 5.9 and A3 (about 0.82)
  X      <- model.matrix(m)
  y      <- model.response(model.frame(m))
  fe_col <- c(1, grep("^factor\\(", colnames(X)))
  r2w    <- 1 - sum(resid(m)^2) / sum(lm.fit(X[, fe_col, drop = FALSE], y)$residuals^2)
  cy_col <- c(1, grep("^factor\\(city_code\\)[^:]*$", colnames(X)))
  r2c    <- 1 - sum(resid(m)^2) / sum(lm.fit(X[, cy_col, drop = FALSE], y)$residuals^2)
  # se_k / p_k: same clustered SE, but the small-sample factor counts ONLY the slope
  # coefficients (not the year dummies). This is the rule behind the SEs printed in
  # Tables 5.10, 5.11, 5.12 and A4 Panels B-C (e.g. 0.080 instead of 0.082 for the full sample).
  k_slope <- sum(!is.na(coef(m)) & !grepl("^\\(Intercept\\)|^factor\\(", names(coef(m))))
  K_used  <- if ("factor(city_code)" %in% attr(terms(m), "term.labels")) m$rank - (vc$G - 1) else m$rank
  se_k    <- se * sqrt((vc$N - K_used) / (vc$N - k_slope))
  p_k     <- 2 * pt(-abs(b / se_k), df = vc$G - 1)
  data.frame(model = label, term = terms, est = unname(b), se = unname(se),
             t = unname(t), p = unname(p), se_k = unname(se_k), p_k = unname(p_k),
             n = vc$N, G = vc$G, r2_within = r2w, r2_city = r2c, row.names = NULL)
}

# =====================================================================
# (3.12) Baseline TWFE
#   HQD_it = b*Dg_it + g'X_it + mu_i + lambda_t + e_it
# =====================================================================
m_312 <- lm(HQD_ent ~ Dg + Con1 + Con2 + Con3 + Con4 + Con5 +
              factor(city_code) + factor(year), data = dat)
res$base <- report(m_312, "3.12 baseline TWFE", c("Dg", CTRL))

# =====================================================================
# (3.13) Subsystem models, entered one at a time
#   HQD_it = b_s*U_s,it + g'X_it + mu_i + lambda_t + e_it
# =====================================================================
m_313T <- lm(HQD_ent ~ UT_u + Con1 + Con2 + Con3 + Con4 + Con5 +
               factor(city_code) + factor(year), data = dat)
m_313P <- lm(HQD_ent ~ UP_u + Con1 + Con2 + Con3 + Con4 + Con5 +
               factor(city_code) + factor(year), data = dat)
m_313I <- lm(HQD_ent ~ UI_u + Con1 + Con2 + Con3 + Con4 + Con5 +
               factor(city_code) + factor(year), data = dat)
res$sub <- rbind(report(m_313T, "3.13 subsystem T", "UT_u"),
                 report(m_313P, "3.13 subsystem P", "UP_u"),
                 report(m_313I, "3.13 subsystem I", "UI_u"))

# =====================================================================
# (3.14) Post-2020 interaction (Post2020 level absorbed by year dummies)
#   HQD_it = b1*Dg + b2*(Dg x Post2020) + g'X + mu_i + lambda_t + e
# =====================================================================
m_314 <- lm(HQD_ent ~ Dg + Dg_post + Con1 + Con2 + Con3 + Con4 + Con5 +
              factor(city_code) + factor(year), data = dat)
res$post <- report(m_314, "3.14 Dg x Post2020", c("Dg", "Dg_post", CTRL))

# post-2020 association b1 + b2 and its clustered SE
vc <- vcov_city(m_314)
b1 <- coef(m_314)[["Dg"]]; b2 <- coef(m_314)[["Dg_post"]]
se_sum <- sqrt(vc$V["Dg", "Dg"] + vc$V["Dg_post", "Dg_post"] + 2 * vc$V["Dg", "Dg_post"])
res$post_sum <- data.frame(model = "3.14 post-2020 association b1+b2", term = "Dg+Dg_post",
                           est = b1 + b2, se = se_sum, t = (b1 + b2) / se_sum,
                           p = 2 * pt(-abs((b1 + b2) / se_sum), vc$G - 1), se_k = NA, p_k = NA,
                           n = vc$N, G = vc$G, r2_within = res$post$r2_within[1],
                           r2_city = res$post$r2_city[1])

# ---- wild cluster bootstrap (Sec 3.4.3(4)) ----------------------------
# The bootstrap uses the twfe()/boot_row() engine from script 7 (not a regression
# formula of its own); it reproduces the same coefficients as m_312 and m_314.
m_base <- twfe(dat, "HQD_ent", c("Dg", CONTROLS))
m_post <- twfe(dat, "HQD_ent", c("Dg", "Dg_post", CONTROLS))
cat("\nWild cluster bootstrap (B =", CFG$boot_B, ", Webb weights, seed", CFG$seed, ")...\n")
m_base_eq  <- twfe(dat, "HQD_eq",  c("Dg", CONTROLS))               # Table 5.8 row (2)
m_post_pca <- twfe(dat, "HQD_pca", c("Dg", "Dg_post", CONTROLS))    # Table 5.8 row (4)
res$boot <- rbind(boot_row(m_base,     "Dg",      "(1) 3.12 baseline, entropy: Dg"),
                  boot_row(m_base_eq,  "Dg",      "(2) 3.12 baseline, equal weight: Dg"),
                  boot_row(m_post,     "Dg_post", "(3) 3.14 entropy: Dg x Post2020"),
                  boot_row(m_post_pca, "Dg_post", "(4) 3.14 PCA: Dg x Post2020"))
save_csv(res$boot, "results_main_inference.csv")

# =====================================================================
# Alternative HQD outcomes in (3.12) and (3.14)
# =====================================================================
m_eq_312  <- lm(HQD_eq ~ Dg + Con1 + Con2 + Con3 + Con4 + Con5 +
                  factor(city_code) + factor(year), data = dat)
m_eq_314  <- lm(HQD_eq ~ Dg + Dg_post + Con1 + Con2 + Con3 + Con4 + Con5 +
                  factor(city_code) + factor(year), data = dat)
m_pca_312 <- lm(HQD_pca ~ Dg + Con1 + Con2 + Con3 + Con4 + Con5 +
                  factor(city_code) + factor(year), data = dat)
m_pca_314 <- lm(HQD_pca ~ Dg + Dg_post + Con1 + Con2 + Con3 + Con4 + Con5 +
                  factor(city_code) + factor(year), data = dat)
res$alt <- rbind(report(m_eq_312,  "baseline, HQD_eq",  "Dg"),
                 report(m_eq_314,  "post2020, HQD_eq",  c("Dg", "Dg_post")),
                 report(m_pca_312, "baseline, HQD_pca", "Dg"),
                 report(m_pca_314, "post2020, HQD_pca", c("Dg", "Dg_post")))

# =====================================================================
# Leave-one-city-out: (3.12) re-estimated without each city / the two cores
# =====================================================================
nm16 <- CITY_NAME[as.character(CITY16)]
nm16[is.na(nm16)] <- as.character(CITY16)[is.na(nm16)]
drops <- c(as.list(CITY16), list(510100, 500000, c(500000, 510100)))
names(drops) <- c(nm16, "Excl. Chengdu", "Excl. Chongqing", "Excl. both cores")
res$loco <- NULL
for (yv in c("HQD_ent", "HQD_eq")) {                  # Table 5.10 / A4B: entropy and equal-weight columns
  for (nm in names(drops)) {
    drop_i <- drops[[nm]]
    m <- lm(as.formula(paste(yv, "~ Dg + Con1 + Con2 + Con3 + Con4 + Con5 +",
                             "factor(city_code) + factor(year)")),
            data = dat[!dat$city_code %in% drop_i, ])
    res$loco <- rbind(res$loco, report(m, paste(yv, "drop", nm), "Dg"))
  }
}

# =====================================================================
# Leave-one-dimension-out HQD outcomes in (3.12)
# =====================================================================
res$lodo <- NULL
for (yv in grep("^HQD_lodo_", names(dat), value = TRUE)) {
  m <- lm(as.formula(paste(yv, "~ Dg + Con1 + Con2 + Con3 + Con4 + Con5 +",
                           "factor(city_code) + factor(year)")), data = dat)
  res$lodo <- rbind(res$lodo, report(m, paste("LODO", yv), "Dg"))
}

# =====================================================================
# (3.15) City-specific linear trends
#   HQD_it = b*Dg + g'X + mu_i + lambda_t + theta_i*t + e
#   t = 1 in 2011, 2 in 2012, ..., 13 in 2023.
#   factor(city_code):t gives every city its own straight-line trend theta_i*t:
#   R creates one column per city that equals t in that city's rows and 0 in
#   all other rows. One of these 16 columns is shown as NA because the year
#   dummies already contain the common trend; this does not affect Dg.
# =====================================================================
dat$t <- dat$year - 2010
res$trend <- NULL
for (yv in c("HQD_ent", "HQD_eq")) {                  # Table A4 Panel A: both HQD measures
  m_315 <- lm(as.formula(paste(yv, "~ Dg + Con1 + Con2 + Con3 + Con4 + Con5 +",
                               "factor(city_code) + factor(year) + factor(city_code):t")), data = dat)
  res$trend <- rbind(res$trend, report(m_315, paste("3.15 city-specific trends,", yv), "Dg"))
}

# =====================================================================
# (3.16) First differences (differencing removes mu_i: no city dummies)
#   dHQD_it = b*dDg_it + g'dX_it + lambda_t + de_it
# =====================================================================
dd <- dat[order(dat$city_code, dat$year), ]
for (v in c("HQD_ent", "HQD_eq", "Dg", CTRL))
  dd[[paste0("d_", v)]] <- ave(dd[[v]], dd$city_code, FUN = function(z) c(NA, diff(z)))
dd <- dd[dd$year > min(dd$year), ]
res$fd <- NULL
for (yv in c("HQD_ent", "HQD_eq")) {
  m_316 <- lm(as.formula(paste0("d_", yv, " ~ d_Dg + d_Con1 + d_Con2 + d_Con3 + d_Con4 + d_Con5 + factor(year)")),
              data = dd)
  res$fd <- rbind(res$fd, report(m_316, paste("3.16 first differences,", yv), "d_Dg"))
}

# =====================================================================
# (3.17)-(3.20) Strength, CV, Balance, Gap
# =====================================================================
U3 <- as.matrix(dat[, c("UT_u", "UP_u", "UI_u")])
dat$Strength <- rowMeans(U3)                                                   # (3.17)
dat$CV       <- sqrt(rowMeans((U3 - rowMeans(U3))^2)) / dat$Strength           # (3.18) divisor 3
dat$Balance  <- 1 - (dat$CV - min(dat$CV)) / (max(dat$CV) - min(dat$CV))       # (3.19)
dat$Gap <- NA
for (i in 1:nrow(dat)) {
  u <- c(dat$UT_u[i], dat$UP_u[i], dat$UI_u[i])
  dat$Gap[i] <- max(u) - min(u)
}

# (3.21) and (3.22), plus strength-only and balance-only
res$sb <- NULL
for (yv in c("HQD_ent", "HQD_eq")) {                  # Table 5.12 (entropy) and A4 Panel C (equal weight)
  fe <- "+ Con1 + Con2 + Con3 + Con4 + Con5 + factor(city_code) + factor(year)"
  m_S   <- lm(as.formula(paste(yv, "~ Strength", fe)), data = dat)
  m_B   <- lm(as.formula(paste(yv, "~ Balance", fe)), data = dat)
  m_321 <- lm(as.formula(paste(yv, "~ Strength + Balance", fe)), data = dat)
  m_322 <- lm(as.formula(paste(yv, "~ Strength + Gap", fe)), data = dat)
  res$sb <- rbind(res$sb,
                  report(m_S,   paste(yv, "strength only"),         "Strength"),
                  report(m_B,   paste(yv, "balance only"),          "Balance"),
                  report(m_321, paste(yv, "3.21 strength+balance"), c("Strength", "Balance")),
                  report(m_322, paste(yv, "3.22 strength+gap"),     c("Strength", "Gap")))
}

# =====================================================================
# (3.23)-(3.26) Industrial-restructuring conditioning
# =====================================================================
# change variables: the thesis uses a FROZEN file of annual changes that already
# contains the 2010 -> 2011 change, so the lagged models have N = 192 (2012-2023).
# Any file in C:/R/spatial whose name ends in "conditioning_variables_final.csv" is used.
cond_file <- list.files(pattern = "conditioning_variables_final\\.csv$")[1]
if (!is.na(cond_file)) {
  cv <- read.csv(cond_file)
  cv <- cv[, c("city_code", "year", "ind_upgrade_change", "ter_share_change", "sec_ter_ratio_change")]
  names(cv) <- c("city_code", "year", "dI1", "dI2", "dRatio")   # dI1 = upgrading, dI2 = tertiary share, dRatio = sec/ter ratio
  dat <- dat[, setdiff(names(dat), c("dI1", "dI2", "dRatio"))]
  dat <- merge(dat, cv, by = c("city_code", "year"))
  dat <- dat[order(dat$city_code, dat$year), ]
  stopifnot(nrow(dat) == 208, !anyNA(dat[, c("dI1", "dI2", "dRatio")]))
  cat("\nConditioning variables read from", cond_file, "(N = 192 in Tables 5.13-5.14)\n")
} else {
  message("conditioning_variables_final.csv not found: changes computed from I1, I2, Con2/I2 (N = 176, not 192)")
}
dat$ratio <- dat$Con2 / dat$I2                                   # secondary / tertiary (used only if no file)
if (!"dI1"    %in% names(dat)) dat$dI1    <- ave(dat$I1,    dat$city_code, FUN = function(z) c(NA, diff(z)))
if (!"dI2"    %in% names(dat)) dat$dI2    <- ave(dat$I2,    dat$city_code, FUN = function(z) c(NA, diff(z)))
if (!"dRatio" %in% names(dat)) dat$dRatio <- ave(dat$ratio, dat$city_code, FUN = function(z) c(NA, diff(z)))
RESTR <- c("dI1", "dI2", "dRatio")

key <- paste(dat$city_code, dat$year)
for (M in RESTR) {
  # within-year z-score with population SD
  z <- rep(NA_real_, nrow(dat))
  for (yr in unique(dat$year)) {
    i <- dat$year == yr
    x <- dat[[M]][i]
    if (!anyNA(x)) z[i] <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  }
  dat[[paste0("Z_", M)]]  <- z
  dat[[paste0("LZ_", M)]] <- z[match(paste(dat$city_code, dat$year - 1), key)]   # (3.23)
}

# (3.27)-(3.29) transport-population index (industry excluded)
dat$T_TP <- (dat$UT_u + dat$UP_u) / 2                                               # (3.27)
dat$C_TP <- ifelse(dat$UT_u + dat$UP_u == 0, 0,
                   2 * sqrt(dat$UT_u * dat$UP_u) / (dat$UT_u + dat$UP_u))           # (3.28)
dat$TPg  <- sqrt(dat$C_TP * dat$T_TP)                                               # (3.29)

res$cond <- NULL; res$slopes <- NULL
for (M in RESTR) {
  lz <- paste0("LZ_", M)
  d  <- dat[!is.na(dat[[lz]]), ]                    # common sample for both models
  # (3.24) and (3.25): centre within the current year
  d$Dc  <- d$Dg   - ave(d$Dg,   d$year)             # (3.24) D^c
  d$TPc <- d$TPg  - ave(d$TPg,  d$year)
  d$LZc <- d[[lz]] - ave(d[[lz]], d$year)           # (3.25) LZ(M)^c
  d$LZ     <- d[[lz]]
  d$INT_D  <- d$Dc  * d$LZc
  d$INT_TP <- d$TPc * d$LZc

  for (yv in c("HQD_ent", "HQD_eq")) {               # Tables 5.13 / 5.14: Panel A entropy, Panel B equal weight
    fe <- "+ Con1 + Con2 + Con3 + Con4 + Con5 + factor(city_code) + factor(year)"
    # (3.25a) HQD = b1 Dg + b2 LZ(M) + b3 [D^c x LZ(M)^c] + g'X + mu_i + lambda_t + e
    m_full <- lm(as.formula(paste(yv, "~ Dg + LZ + INT_D", fe)), data = d)
    # (3.30)  HQD = t1 TP + t2 LZ(M) + t3 [TP^c x LZ(M)^c] + g'X + mu_i + lambda_t + e
    m_tp   <- lm(as.formula(paste(yv, "~ TPg + LZ + INT_TP", fe)), data = d)
    r_full <- report(m_full, paste(yv, "3.25 full TPI x", M), c("Dg", "LZ", "INT_D"))
    r_tp   <- report(m_tp,   paste(yv, "3.30 TP-only x",  M), c("TPg", "LZ", "INT_TP"))
    r_full$term <- c("Dg", lz, paste0("INT_Dg_", M))
    r_tp$term   <- c("TPg", lz, paste0("INT_TPg_", M))
    res$cond <- rbind(res$cond, r_full, r_tp)

    # (3.26) conditional slope b1 + b3 * LZ(M)^c at -1, 0, +1 SD
    vc <- vcov_city(m_full)
    for (k in c(-1, 0, 1)) {
      est <- coef(m_full)[["Dg"]] + k * coef(m_full)[["INT_D"]]
      se  <- sqrt(vc$V["Dg", "Dg"] + k^2 * vc$V["INT_D", "INT_D"] + 2 * k * vc$V["Dg", "INT_D"])
      res$slopes <- rbind(res$slopes, data.frame(outcome = yv, restructuring = M, LZc = k, slope = est,
                                                 se = se, p = 2 * pt(-abs(est / se), vc$G - 1), n = vc$N))
    }
  }
}

# =====================================================================
# (3.31)-(3.37) spatial and network weight matrices (base R, no packages)
# ---------------------------------------------------------------------
# Inputs (in C:/R/spatial). If a file is missing, that matrix is skipped;
# nothing is simulated.
#   city_centroids_3857.csv : city_code, x, y   (polygon centroids, EPSG:3857, metres)
#   queen_neighbours.csv    : city_code_i, city_code_j   (one row per neighbouring pair)
#   hsr_2015.csv, hsr_2019.csv, hsr_2021.csv, hsr_2023.csv :
#       16 x 16 shortest HSR minutes; first column = city_code, header = city codes;
#       leave a cell BLANK when there is no service (never 0)
# =====================================================================
city_order <- sort(unique(dat$city_code)); n_city <- length(city_order)
W_ALL <- list()

# (3.31)-(3.33): w_ii = 0, row-standardise non-empty rows, check row sums = 1
row_std <- function(W0) {
  diag(W0) <- 0                                                  # (3.31)
  rs <- rowSums(W0); W <- W0
  W[rs > 0, ] <- W0[rs > 0, , drop = FALSE] / rs[rs > 0]         # (3.32)
  stopifnot(all(abs(rowSums(W)[rs > 0] - 1) < 1e-10))            # (3.33)
  if (any(rs == 0)) message("  note: ", sum(rs == 0), " city row(s) have no link (left empty)")
  W
}

if (file.exists("city_centroids_3857.csv")) {
  cen_xy <- read.csv("city_centroids_3857.csv")
  cen_xy <- cen_xy[match(city_order, cen_xy$city_code), ]
  xy     <- as.matrix(cen_xy[, c("x", "y")])
  dist_m <- as.matrix(dist(xy))                                  # d_ij in metres

  # (3.35) four nearest neighbours (not symmetrised, as in libpysal)
  W0 <- matrix(0, n_city, n_city)
  for (i in seq_len(n_city)) { di <- dist_m[i, ]; di[i] <- Inf; W0[i, order(di)[1:4]] <- 1 }
  W_ALL$knn4 <- row_std(W0)

  # (3.36) inverse distance: 1/d_ij for i != j
  W0 <- 1 / dist_m; diag(W0) <- 0
  W_ALL$idw <- row_std(W0)
} else message("city_centroids_3857.csv not found: (3.35) and (3.36) skipped")

if (file.exists("queen_neighbours.csv")) {
  qp <- read.csv("queen_neighbours.csv")
  W0 <- matrix(0, n_city, n_city)
  ii <- match(qp$city_code_i, city_order); jj <- match(qp$city_code_j, city_order)
  W0[cbind(ii, jj)] <- 1; W0[cbind(jj, ii)] <- 1                 # (3.34) shared boundary or vertex
  W_ALL$queen <- row_std(W0)
} else message("queen_neighbours.csv not found: (3.34) skipped")

# (3.37) HSR inverse travel time: 1/tau_ij,s if a valid time is observed, else 0
for (s in c("2015", "2019", "2021", "2023")) {
  f <- paste0("hsr_", s, ".csv")
  if (!file.exists(f)) { message(f, " not found: (3.37) for ", s, " skipped"); next }
  h   <- read.csv(f, check.names = FALSE, na.strings = c("", "NA", "-"))
  tau <- as.matrix(h[, -1]); rownames(tau) <- h[[1]]
  tau <- tau[as.character(city_order), as.character(city_order)]
  storage.mode(tau) <- "numeric"
  W0  <- ifelse(!is.na(tau) & tau > 0, 1 / tau, 0)
  W_ALL[[paste0("hsr_", s)]] <- row_std(W0)
}

for (nm in names(W_ALL)) {
  save_csv(data.frame(city_code = city_order, W_ALL[[nm]], check.names = FALSE),
           paste0("W_", nm, ".csv"))
  dimnames(W_ALL[[nm]]) <- list(city_order, city_order)
}
saveRDS(W_ALL, outp("weight_matrices.rds"))
cat("\nWeight matrices built:", if (length(W_ALL)) paste(names(W_ALL), collapse = ", ") else "none", "\n")

# =====================================================================
# PART C  PRINT AND SAVE ALL RESULTS
# =====================================================================
labels <- c(
  base     = "(3.12) Baseline TWFE",
  sub      = "(3.13) Subsystem models (T, P, I)",
  post     = "(3.14) Dg x Post2020",
  post_sum = "(3.14) Post-2020 slope b1 + b2",
  boot     = "Wild cluster bootstrap for (3.12) and (3.14)",
  alt      = "Alternative HQD outcomes: (3.12) and (3.14)",
  loco     = "Leave-one-city-out: (3.12)",
  lodo     = "Leave-one-dimension-out: (3.12)",
  trend    = "(3.15) City-specific linear trends",
  fd       = "(3.16) First differences",
  sb       = "(3.21)-(3.22) Strength-balance decomposition",
  cond     = "(3.25a) Industrial conditioning and (3.30) TP-only diagnostic",
  slopes   = "(3.26) Conditional slopes at LZc = -1, 0, +1")

for (k in names(labels)) if (!is.null(res[[k]])) {
  cat("\n==========", labels[[k]], "==========\n")
  print(res[[k]], digits = 4, row.names = FALSE)
  save_csv(res[[k]], paste0("results_", k, ".csv"))
}

# all regression tables in one file (bootstrap and slopes have other columns)
all_res <- do.call(rbind, lapply(setdiff(names(labels), c("boot", "slopes")), function(k) {
  x <- res[[k]]; if (is.null(x)) return(NULL); x$block <- k; x }))
save_csv(all_res, "results_all_models.csv")

# Constructed variables: (3.17)-(3.20), (3.23), (3.27)-(3.29)
vars <- intersect(c("Strength", "CV", "Balance", "Gap", "T_TP", "C_TP", "TPg",
                    paste0("LZ_", RESTR)), names(dat))
desc <- data.frame(variable = vars,
                   n    = sapply(dat[vars], function(x) sum(!is.na(x))),
                   mean = sapply(dat[vars], mean, na.rm = TRUE),
                   sd   = sapply(dat[vars], sd,   na.rm = TRUE),
                   min  = sapply(dat[vars], min,  na.rm = TRUE),
                   max  = sapply(dat[vars], max,  na.rm = TRUE))
cat("\n========== Constructed variables (3.17)-(3.20), (3.23), (3.27)-(3.29) ==========\n")
print(desc, digits = 4, row.names = FALSE)
save_csv(desc, "results_constructed_variables.csv")
save_csv(dat,  "analysis_panel_all_variables.csv")
saveRDS(res, outp(tag("results_R3.rds"))); saveRDS(dat, outp(tag("analysis_panel.rds")))

# =====================================================================
# TABLE 5.1  Summary statistics of variables (N, Mean, SD, Min, Max)
#   SD is the sample standard deviation (divides by n - 1), the same as
#   Stata's "summarize". The last columns compare each value with the thesis.
# =====================================================================
T51_VARS  <- c("HQD_ent", "HQD_eq", "HQD_pca", "Dg", "Con1", "Con2", "Con3", "Con4", "Con5")
T51_LABEL <- c("hqd_entropy", "hqd_equal", "hqd_pca", "Dg", "Con1", "Con2", "Con3", "Con4", "Con5")

# Values printed in Table 5.1 of the thesis (used only for the comparison)
T51_THESIS <- data.frame(
  mean = c(0.1874, 0.3762, 0.3152, 0.5548, 10.6438, 0.4260, 7.8649, 0.1939, 0.0325),
  sd   = c(0.1044, 0.0875, 0.2042, 0.0999,  0.4057, 0.0766, 0.6028, 0.0660, 0.0098),
  min  = c(0.0387, 0.1939, 0.0000, 0.3646,  9.7000, 0.2800, 6.5300, 0.1100, 0.0200),
  max  = c(0.5649, 0.6225, 1.0000, 0.8629, 11.5500, 0.6100, 9.4400, 0.6800, 0.0900))

tab51 <- data.frame(Variable = T51_LABEL, N = NA, Mean = NA, SD = NA, Min = NA, Max = NA)

for (k in 1:length(T51_VARS)) {
  x <- dat[[T51_VARS[k]]]                # one column of the panel
  if (is.null(x)) next                   # variable not in the data: leave NA
  x <- x[!is.na(x)]
  tab51$N[k]    <- length(x)
  tab51$Mean[k] <- mean(x)
  tab51$SD[k]   <- sd(x)                 # sample SD (n - 1)
  tab51$Min[k]  <- min(x)
  tab51$Max[k]  <- max(x)
}

# Comparison with the thesis: same = all four statistics agree to 4 decimals
tab51$Mean_thesis <- T51_THESIS$mean
tab51$SD_thesis   <- T51_THESIS$sd
tol <- 0.00006                           # allows for rounding in the thesis
tab51$Same_as_thesis <- ifelse(abs(tab51$Mean - T51_THESIS$mean) < tol &
                               abs(tab51$SD   - T51_THESIS$sd)   < tol &
                               abs(tab51$Min  - T51_THESIS$min)  < tol &
                               abs(tab51$Max  - T51_THESIS$max)  < tol, "yes", "NO")

cat("\n========== Table 5.1 Summary statistics of variables ==========\n")
if (DEMO_EXTRA) cat("NOTE: hqd_controls.xlsx not found, so HQD and Con1-Con5 are SIMULATED.\n",
                    "     Only the Dg row can be compared with the thesis.\n")
print(data.frame(tab51[, 1:2], round(tab51[, 3:8], 4), Same_as_thesis = tab51$Same_as_thesis),
      row.names = FALSE)
save_csv(tab51, "table_5_1_summary_statistics.csv")

# Weight matrices (3.31)-(3.37)
cat("\n========== (3.31)-(3.37) Weight matrices ==========\n")
if (length(W_ALL)) {
  for (nm in names(W_ALL)) { cat("\n--", nm, "--\n"); print(round(W_ALL[[nm]], 3)) }
} else cat("None built: the input files are not in the folder (see messages above).\n")

cat("\nAll tables saved in:", normalizePath(dirname(outp("x"))), "\n")
