# =====================================================================
# R7  Reproducibility and Computational Implementation   (Section 3.4.7)
# ---------------------------------------------------------------------
# Numbered equations: none. Section 3.4.7 prescribes a workflow, so this
# file implements it: configuration, fixed seeds, panel-key validation,
# no-imputation rule, input hashing, session log, and a master runner.
#
# Every other script starts with  source("R7_reproducibility.R").
# To run the whole pipeline:      Sys.setenv(TPI_RUN_ALL="1"); source("R7_reproducibility.R")
# =====================================================================
setwd("C:/R/spatial")
suppressPackageStartupMessages({ library(readxl) })
options(scipen = 999)                      # print 500000, not 5e+05 (city labels)

# ---- 1. Configuration ------------------------------------------------
# panel_data_final.xlsx holds BOTH the TPI indicators and the HQD indicators + controls
CFG <- list(
  data_file   = "panel_data_final.xlsx",   # TPI indicators (T1 ... I2)
  extra_file  = "panel_data_final.xlsx",   # HQD indicators (Inno1 ... Sha3) + controls (Con1 ... Con5)
  shape_file  = "ccec_cities.shp",         # OPTIONAL: 16 city polygons, field `city_code`
  hsr_file    = "hsr_travel_time.xlsx",    # OPTIONAL: sheets "2015","2019","2021","2023"
  out_dir     = "output",
  years       = 2011:2023,
  seed        = 20260730,                  # Section 3.4.3(4)
  boot_B      = as.integer(Sys.getenv("TPI_BOOT_B", "9999")),   # wild cluster bootstrap reps
  ci_grid     = as.integer(Sys.getenv("TPI_CI_GRID", "41")),    # grid points for CI inversion
  perm_B      = 999L                       # Moran / LISA permutations
)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# The 16 city units used in the thesis (Chongqing + 15 Sichuan prefecture-level cities).
# The supplied workbook also contains Guangyuan (510800) and Bazhong (511900); they are
# NOT part of the thesis sample (208 obs = 16 x 13) and are dropped here.
CITY16 <- c(500000, 510100, 510300, 510500, 510600, 510700, 510900, 511000,
            511100, 511300, 511400, 511500, 511600, 511700, 511800, 512000)
CITY_NAME <- c(`500000`="Chongqing", `510100`="Chengdu", `510300`="Zigong", `510500`="Luzhou",
               `510600`="Deyang", `510700`="Mianyang", `510900`="Suining", `511000`="Neijiang",
               `511100`="Leshan", `511300`="Nanchong", `511400`="Meishan", `511500`="Yibin",
               `511600`="Guang'an", `511700`="Dazhou", `511800`="Ya'an", `512000`="Ziyang")

# TPI indicator system (Table 3.2); all indicators are positive.
TPI_GROUPS <- list(T = c("T1", "T2"), P = c("lnP1", "P2", "P3"), I = c("I1", "I2"))

set.seed(CFG$seed)
outp <- function(...) file.path(CFG$out_dir, ...)

# ---- 2. Data loading and validation (no interpolation / imputation) ---
validate_panel <- function(df, n_city = 16, years = CFG$years) {
  need <- unlist(TPI_GROUPS)
  miss <- setdiff(c("city_code", "year", need), names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  for (v in need) if (!is.numeric(df[[v]])) stop("Non-numeric TPI indicator: ", v)
  if (anyNA(df[, need])) stop("Missing TPI values found; scripts stop rather than impute (Sec. 3.2.2).")
  if (anyDuplicated(df[, c("city_code", "year")])) stop("Panel key (city_code, year) is not unique.")
  if (length(unique(df$city_code)) != n_city) stop("Expected ", n_city, " cities.")
  if (nrow(df) != n_city * length(years)) stop("Expected ", n_city * length(years), " observations.")
  ok <- tapply(df$year, df$city_code, function(y) identical(sort(y), as.numeric(years)) || identical(sort(y), years))
  if (!all(ok)) stop("Years are not consecutive 2011-2023 within every city.")
  invisible(TRUE)
}

load_tpi_panel <- function() {
  d <- as.data.frame(read_excel(CFG$data_file))
  names(d)[names(d) == "Ln(P1)"] <- "lnP1"
  d <- d[d$city_code %in% CITY16, ]
  d$city_code <- as.numeric(d$city_code); d$year <- as.numeric(d$year)
  d$city <- unname(CITY_NAME[as.character(d$city_code)])
  d <- d[order(d$city_code, d$year), ]
  rownames(d) <- NULL
  validate_panel(d)
  d
}

# ---- 3. Optional inputs (HQD indicators, controls, geometry, HSR) -----
# hqd_controls.xlsx must contain: city_code, year,
#   Inno1 Inno2 Coor1 Coor2 Green1 Green2 Open1 Sha1 Sha2 Sha3   (raw values; Sha1 in levels - ln is applied here)
#   Con1 Con2 Con3 Con4 Con5                                       (Table 3.3)
# optional: dI1 dI2 dRatio  (frozen annual-change variables; if absent they are computed
#   from I1, I2 and Con2/I2, which needs 2011 changes from 2010 data to reach 192 obs)
HQD_IND <- data.frame(
  code = c("Inno1","Inno2","Coor1","Coor2","Green1","Green2","Open1","Sha1","Sha2","Sha3"),
  dim  = c("Innovation","Innovation","Coordination","Coordination","Green","Green","Openness","Sharing","Sharing","Sharing"),
  dir  = c("+","+","-","-","-","+","+","+","+","+"), stringsAsFactors = FALSE)
CONTROLS <- paste0("Con", 1:5)

# City centroids (approx. lon/lat) used ONLY to build demonstration geometry.
.demo_xy <- data.frame(
  city_code = CITY16,
  lon = c(106.55,104.07,104.78,105.44,104.40,104.68,105.59,105.06,103.77,106.11,103.85,104.64,106.63,107.47,103.01,104.63),
  lat = c(29.56,30.67,29.34,28.87,31.13,31.47,30.53,29.58,29.55,30.84,30.08,28.77,30.46,31.21,29.98,30.13))

make_demo_extra <- function(panel) {
  set.seed(CFG$seed + 1)
  n <- nrow(panel); cf <- rnorm(length(CITY16))[match(panel$city_code, CITY16)]
  tt <- panel$year - 2011
  f <- function(base, sd1 = 1) base + 0.3 * cf + 0.05 * tt + rnorm(n, 0, sd1) * 0.3
  e <- data.frame(city_code = panel$city_code, year = panel$year)
  e$Inno1 <- exp(f(0)); e$Inno2 <- exp(f(-1)); e$Coor1 <- exp(f(1) * 0.3) + 1
  e$Coor2 <- exp(f(-2)); e$Green1 <- exp(f(0)); e$Green2 <- 35 + 3 * f(0)
  e$Open1 <- exp(f(-2)); e$Sha1 <- exp(f(8)); e$Sha2 <- 5 + f(0); e$Sha3 <- 90 + 3 * f(0)
  e$Con1 <- f(10); e$Con2 <- 45 + 3 * f(0); e$Con3 <- f(7); e$Con4 <- exp(f(-2)); e$Con5 <- exp(f(-4))
  e
}

DEMO_EXTRA <- !file.exists(CFG$extra_file)
load_extra <- function(panel) {
  if (DEMO_EXTRA) {
    message("*** DEMO MODE: '", CFG$extra_file, "' not found. HQD indicators and controls are SIMULATED. ",
            "Any output file starting with DEMO_ is a code test only and must not be reported. ***")
    return(make_demo_extra(panel))
  }
  e <- as.data.frame(read_excel(CFG$extra_file))
  e$city_code <- as.numeric(e$city_code); e$year <- as.numeric(e$year)
  e <- e[e$city_code %in% CITY16, ]
  need <- c(HQD_IND$code, CONTROLS)
  if (length(setdiff(need, names(e)))) stop("hqd_controls.xlsx lacks: ", paste(setdiff(need, names(e)), collapse = ", "))
  if (anyNA(e[, need])) stop("Missing HQD/control values: not imputed (Sec. 3.3.1).")
  e
}
tag <- function(name) paste0(if (DEMO_EXTRA) "DEMO_" else "", name)

# ---- 3b. Shared estimation engine (used by R3 and R4) -------------------
# ======================= estimation engine ===============================
# OLS with dummy-variable FE and city-clustered covariance.
design <- function(d, xvars, city_fe = TRUE, year_fe = TRUE, city_trend = FALSE) {
  parts <- c(xvars, if (city_fe) "factor(city_code)", if (year_fe) "factor(year)",
             if (city_trend) "factor(city_code):year")
  model.matrix(as.formula(paste("~", paste(parts, collapse = "+"))), data = d)
}
fit_cl <- function(y, X, cl) {
  q <- qr(X); r <- q$rank; keep <- sort(q$pivot[seq_len(r)]); Xr <- X[, keep, drop = FALSE]
  XtXi <- chol2inv(chol(crossprod(Xr))); dimnames(XtXi) <- list(colnames(Xr), colnames(Xr)); b <- drop(XtXi %*% crossprod(Xr, y)); e <- drop(y - Xr %*% b)
  n <- length(y); G <- length(unique(cl))
  k_eff <- r - sum(grepl("^factor\\(city_code\\)[^:]*$", colnames(Xr)))   # city dummies are nested in clusters: not counted in K
  adj <- G / (G - 1) * (n - 1) / (n - k_eff)
  S <- rowsum(Xr * e, cl); V <- adj * XtXi %*% crossprod(S) %*% XtXi
  list(b = setNames(b, colnames(Xr)), V = V, e = e, n = n, G = G, k = r, adj = adj, X = Xr, y = y)
}
twfe <- function(d, y, xvars, ..., cluster = "city_code") {
  d <- d[complete.cases(d[, c(y, xvars)]), ]
  X <- design(d, xvars, ...); f <- fit_cl(d[[y]], X, d[[cluster]])
  f$data <- d; f$xvars <- xvars; f$yname <- y
  # within R2 (variation after removing FE); uses residualised y and fitted slopes
  fe_only <- if (ncol(X) > length(xvars) + 1) X[, setdiff(colnames(X), xvars), drop = FALSE] else X[, 1, drop = FALSE]
  y_res <- resid(lm.fit(fe_only, d[[y]])); f$r2_within <- 1 - sum(f$e^2) / sum(y_res^2)
  f
}
coef_tab <- function(f, label, terms = f$xvars) {
  terms <- intersect(terms, names(f$b)); se <- sqrt(diag(f$V))[terms]; t <- f$b[terms] / se
  # se_k: small-sample factor counts only the slope coefficients (rule used in Tables 5.10-5.15, A4, A5)
  adj_k <- f$G / (f$G - 1) * (f$n - 1) / (f$n - length(f$xvars))
  se_k  <- se * sqrt(adj_k / f$adj)
  data.frame(model = label, term = terms, est = unname(f$b[terms]), se = unname(se), t = unname(t),
             p = 2 * pt(-abs(t), df = f$G - 1), se_k = unname(se_k),
             p_k = unname(2 * pt(-abs(f$b[terms] / se_k), df = f$G - 1)),
             n = f$n, G = f$G, r2_within = f$r2_within, row.names = NULL)
}

# ---- wild cluster bootstrap-t, null imposed, Webb 6-point weights -----------
WEBB <- c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5))
wcr_p <- function(f, test, beta0 = 0, B = CFG$boot_B, seed = CFG$seed) {
  X <- f$X; y <- f$y
  cl <- as.integer(factor(f$data$city_code)); G <- max(cl)
  j <- match(test, colnames(X)); x <- X[, j]; Xo <- X[, -j, drop = FALSE]
  xt <- resid(lm.fit(Xo, x)); sxx <- sum(xt^2)                     # FWL residualised regressor
  XtXi <- chol2inv(chol(crossprod(X)))
  beta_hat <- sum(xt * y) / sxx; e <- drop(y - X %*% (XtXi %*% crossprod(X, y)))
  se_hat <- sqrt(f$adj * sum(rowsum(xt * e, cl)^2) / sxx^2); t0 <- (beta_hat - beta0) / se_hat
  y0 <- y - beta0 * x; fit_r <- lm.fit(Xo, y0)                     # restricted model under H0
  base <- beta0 * x + (y0 - resid(fit_r)); u <- resid(fit_r)
  set.seed(seed); Wg <- matrix(sample(WEBB, G * B, replace = TRUE), G, B)
  Ys <- base + u * Wg[cl, ]                                        # n x B bootstrap outcomes
  Es <- Ys - X %*% (XtXi %*% crossprod(X, Ys))
  bs <- drop(crossprod(xt, Ys)) / sxx
  ses <- sqrt(f$adj * colSums(rowsum(xt * Es, cl)^2) / sxx^2)
  tstar <- (bs - beta0) / ses
  list(p = mean(abs(tstar) >= abs(t0) - 1e-12), t0 = t0, beta = beta_hat, se = se_hat)
}
wcr_ci <- function(f, test, level = 0.95, B = CFG$boot_B) {
  est <- f$b[[test]]; se <- sqrt(f$V[test, test]); alpha <- 1 - level; width <- 5
  for (try in 1:4) {                      # widen the search grid until both ends reject
    grid <- seq(est - width * se, est + width * se, length.out = CFG$ci_grid)
    pv <- sapply(grid, function(b0) wcr_p(f, test, beta0 = b0, B = B)$p)
    if (pv[1] <= alpha && pv[length(pv)] <= alpha) break
    width <- width * 2
  }
  inside <- which(pv > alpha)
  if (!length(inside)) return(c(NA_real_, NA_real_))
  cross <- function(i0, i1) grid[i0] + (alpha - pv[i0]) * (grid[i1] - grid[i0]) / (pv[i1] - pv[i0])  # interpolate p = alpha
  lo <- min(inside); hi <- max(inside)
  c(if (lo == 1) NA_real_ else cross(lo, lo - 1), if (hi == length(grid)) NA_real_ else cross(hi, hi + 1))
}
boot_row <- function(f, test, label) {
  b <- wcr_p(f, test); ci95 <- wcr_ci(f, test, 0.95); ci90 <- wcr_ci(f, test, 0.90)
  data.frame(model = label, term = test, est = f$b[[test]], se_cluster = sqrt(f$V[test, test]),
             p_conv = 2 * pt(-abs(f$b[[test]] / sqrt(f$V[test, test])), f$G - 1), p_boot = b$p,
             ci90_lo = ci90[1], ci90_hi = ci90[2], ci95_lo = ci95[1], ci95_hi = ci95[2],
             B = CFG$boot_B, seed = CFG$seed, row.names = NULL)
}


# ---- 4. Run-order helper ----------------------------------------------
# Internal script names -> your file names in C:/R/spatial
SCRIPT_FILE <- c(R1_standardisation_entropy.R        = "1. R_code for entropy.R",
                 R2_coupling_coordination.R          = "2. R_code for coupling_coordination.R",
                 R3_panel_regression_sensitivity.R   = "3. R_code for panel_regression(e).R",
                 R4_spatial_network.R                = "4. R_code for spatial_network(b).R",
                 R5_distribution_dynamics.R          = "5. R_code for distribution_dynamics(b).R",
                 R6_evidence_hierarchy.R             = "6. R_code for evidence_hierarchy.R")
need_output <- function(file, script) {
  if (!file.exists(outp(file))) {
    if (script %in% names(SCRIPT_FILE)) script <- SCRIPT_FILE[[script]]
    message("Running ", script, " to create ", file); source(script, local = FALSE)
  }
}

# ---- 5. Manifest: input hashes + session information -------------------
write_manifest <- function() {
  inputs <- c(CFG$data_file, CFG$extra_file, CFG$shape_file, CFG$hsr_file)
  inputs <- inputs[file.exists(inputs)]
  man <- data.frame(file = inputs, md5 = unname(tools::md5sum(inputs)))
  write.csv(man, outp("manifest_inputs.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()),
               paste("seed:", CFG$seed), paste("boot_B:", CFG$boot_B), paste("perm_B:", CFG$perm_B),
               paste("demo_mode:", DEMO_EXTRA)), outp("session_info.txt"))
  invisible(man)
}

# ---- 6. Master runner ---------------------------------------------------
run_all <- function() {
  scripts <- unname(SCRIPT_FILE)
  for (s in scripts) { message("\n=== ", s, " ==="); source(s, local = FALSE) }
  write_manifest()
  message("\nPipeline finished. Outputs in '", CFG$out_dir, "/'.")
}
if (identical(Sys.getenv("TPI_RUN_ALL"), "1") && !exists(".TPI_RUNNING", envir = globalenv())) {
  assign(".TPI_RUNNING", TRUE, envir = globalenv()); run_all()
}
