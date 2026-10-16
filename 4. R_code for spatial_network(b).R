# =====================================================================
# R4  Spatial and Network Analytical Framework             (Section 3.4.4)
# ---------------------------------------------------------------------
# Equations implemented
#   (3.31) zero diagonal                      (3.32) row standardisation
#   (3.33) row sums = 1 (checked)
#   (3.34) queen contiguity                   (3.35) 4-nearest-neighbour (not symmetrised)
#   (3.36) inverse distance (EPSG:3857 centroids)
#   (3.37) dynamic HSR inverse-time weights   (3.38) snapshot -> year assignment
#   (3.39) global Moran's I     (3.40) S0
#   (3.41) local Moran's I      (3.42) z-score          999 permutations
#   (3.43) contemporaneous exposure  W D     (3.44) lagged exposure
#   (3.45) SLX BASE / LITE / FULL specifications
#   (3.46) within-threshold HSR matrix        (3.47) beyond-threshold HSR matrix
#   multiple testing: Holm and Benjamini-Hochberg within each outcome family
#
# INPUT: real geometry/HSR files if present (ccec_cities.shp, hsr_travel_time.xlsx);
#        otherwise DEMONSTRATION geometry (approximate centroids, Delaunay "queen"
#        stand-in, simulated travel times). Demo outputs are prefixed DEMO_.
# =====================================================================
setwd("C:/R/spatial") 
source("7. R_code for reproducibility.R")
suppressPackageStartupMessages({ library(spdep); library(sf) })
need_output("tpi_index.rds", "R2_coupling_coordination.R")
need_output(tag("analysis_panel.rds"), "R3_panel_regression_sensitivity.R")
dat <- readRDS(outp(tag("analysis_panel.rds")))
# VERSION (b): no simulated geometry any more.
#  * HSR matrices are read from HSR_matrices.xlsx (sheets 2015, 2019, 2021, 2023).
#  * Queen contiguity: from ccec_cities.shp if present; otherwise from the edge list
#    QUEEN_EDGES below (16 city units, shared boundary or vertex). This list reproduces
#    every queen Moran's I in Table A1 (2011-2023) to three decimals.
#  * 4-nearest-neighbour and inverse-distance matrices need polygon centroids, so they
#    are built ONLY when ccec_cities.shp is present; otherwise they are skipped.
CFG$hsr_file <- "HSR_matrices.xlsx"
HAVE_SHAPE <- file.exists(CFG$shape_file)
HAVE_HSR   <- file.exists(CFG$hsr_file)
if (!HAVE_HSR) stop("HSR_matrices.xlsx not found in ", getwd())
if (!HAVE_SHAPE) message("ccec_cities.shp not found: queen matrix taken from QUEEN_EDGES; 4-NN and inverse-distance skipped.")
DEMO_SPATIAL <- FALSE
tag_s <- function(name) paste0(if (DEMO_EXTRA) "DEMO_" else "", name)
QUEEN_EDGES <- rbind(                                   # pairs of city codes that share a boundary
  c(510100,510600), c(510100,511400), c(510100,512000), c(510100,511800), c(510600,510700),
  c(510600,510900), c(510700,510900), c(510700,511300), c(510900,511300), c(510900,511600),
  c(510900,500000), c(510900,512000), c(511300,511600), c(511300,511700), c(511600,511700),
  c(511600,500000), c(511700,500000), c(512000,511000), c(511000,510300), c(511000,510500),
  c(511000,500000), c(510300,510500), c(510300,511500), c(510300,511100), c(510500,511500),
  c(510500,500000), c(511500,511100), c(511100,511400), c(511100,511800), c(511400,511800),
  c(512000,511400), c(511000,511400), c(500000,512000), c(510600,512000), c(510300,511400))
N16 <- length(CITY16); dn <- list(as.character(CITY16), as.character(CITY16))

# ======================= weight-matrix equations ===============================
w_zero_diag <- function(W0) { diag(W0) <- 0; W0 }                                    # (3.31)
row_standardise <- function(W0) {                                                    # (3.32)
  W0 <- w_zero_diag(W0); rs <- rowSums(W0); W <- W0; nz <- rs > 0
  W[nz, ] <- W0[nz, , drop = FALSE] / rs[nz]
  stopifnot(all(abs(rowSums(W)[nz] - 1) < 1e-12))                                    # (3.33)
  W
}
w_queen <- function(nb) { W <- nb2mat(nb, style = "B", zero.policy = TRUE); dimnames(W) <- dn; W }   # (3.34)
w_knn   <- function(xy, k = 4) { nb <- knn2nb(knearneigh(xy, k = k))                 # (3.35) as libpysal: no symmetrisation
                                 W <- nb2mat(nb, style = "B"); dimnames(W) <- dn; W }
w_idw   <- function(xy) { d <- as.matrix(dist(xy)); W <- 1 / d; diag(W) <- 0; dimnames(W) <- dn; W }   # (3.36)
w_hsr_raw <- function(tau) { W <- ifelse(!is.na(tau) & tau > 0, 1 / tau, 0); diag(W) <- 0; dimnames(W) <- dn; W }  # (3.37)
hsr_snapshot <- function(t) {                                                        # (3.38)
  if (t >= 2015 && t <= 2018) "2015" else if (t >= 2019 && t <= 2020) "2019"
  else if (t >= 2021 && t <= 2022) "2021" else if (t == 2023) "2023" else NA_character_
}
w_hsr_within <- function(tau, cc) { W <- ifelse(!is.na(tau) & tau > 0 & tau <= cc, 1 / tau, 0); diag(W) <- 0; dimnames(W) <- dn; W }  # (3.46)
w_hsr_beyond <- function(tau, cc) { W <- ifelse(!is.na(tau) & tau > cc, 1 / tau, 0);            diag(W) <- 0; dimnames(W) <- dn; W }  # (3.47)

# ======================= inputs: geometry and HSR ==============================
# HSR travel times (3.37): blank cells stay NA (= no service, never zero minutes)
hsr_tau <- lapply(c("2015", "2019", "2021", "2023"), function(s) {
  m <- as.data.frame(read_excel(CFG$hsr_file, sheet = s)); ids <- as.numeric(m[[1]]); m <- as.matrix(m[, -1])
  m <- m[match(CITY16, ids), match(CITY16, as.numeric(colnames(m)))]
  storage.mode(m) <- "numeric"; dimnames(m) <- dn; m })
names(hsr_tau) <- c("2015", "2019", "2021", "2023")

xy <- NULL
if (HAVE_SHAPE) {
  g <- st_read(CFG$shape_file, quiet = TRUE); g <- g[match(CITY16, g$city_code), ]
  stopifnot(nrow(g) == N16, !anyNA(g$city_code))
  Wq0 <- w_queen(poly2nb(g, queen = TRUE))
  xy  <- st_coordinates(st_centroid(st_geometry(st_transform(g, 3857))))
} else {
  Wq0 <- matrix(0, N16, N16, dimnames = dn)                  # (3.34) from the edge list
  for (k in seq_len(nrow(QUEEN_EDGES))) {
    i <- match(QUEEN_EDGES[k, 1], CITY16); j <- match(QUEEN_EDGES[k, 2], CITY16)
    Wq0[i, j] <- 1; Wq0[j, i] <- 1
  }
}
Wq <- row_standardise(Wq0)
stopifnot(all(rowSums(Wq) > 0))                              # every unit has >=1 queen neighbour
Wfun <- list(queen = function(t) Wq)
if (!is.null(xy)) {
  Wk <- row_standardise(w_knn(xy, 4)); Wd <- row_standardise(w_idw(xy))
  Wfun$knn4 <- function(t) Wk; Wfun$idw <- function(t) Wd
}
Wfun$hsr <- function(t) { s <- hsr_snapshot(t); if (is.na(s)) NULL else row_standardise(w_hsr_raw(hsr_tau[[s]])) }
saveRDS(list(queen = Wq, hsr_tau = hsr_tau, xy = xy), outp(tag_s("spatial_weights.rds")))

# ======================= Moran's I and LISA ====================================
moran_I <- function(x, W) {                                  # (3.39) with S0 from (3.40)
  z <- x - mean(x); N <- length(x); S0 <- sum(W)
  N / S0 * sum(W * outer(z, z)) / sum(z^2)
}
local_moran <- function(x, W) { z <- (x - mean(x)) / sqrt(mean((x - mean(x))^2)); as.vector(z * (W %*% z)) }  # (3.41),(3.42)
p_sim <- function(obs, sims) { B <- length(sims); larger <- sum(sims >= obs); if (B - larger < larger) larger <- B - larger; (larger + 1) / (B + 1) }
active_W <- function(W0) {                                   # drop cities with no positive weight, re-standardise
  keep <- which(rowSums(W0) > 0); list(keep = keep, W = row_standardise(W0[keep, keep, drop = FALSE])) }

moran_year <- function(x, W, B = CFG$perm_B) {
  I <- moran_I(x, W); sims <- replicate(B, moran_I(sample(x), W))
  li <- local_moran(x, W); N <- length(x); z <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  lp <- sapply(seq_len(N), function(i) {
    nz <- which(W[i, ] > 0); if (!length(nz)) return(NA_real_)
    oth <- setdiff(seq_len(N), i)
    sims_i <- replicate(B, { idx <- sample(oth, length(nz)); z[i] * sum(W[i, nz] * z[idx]) }); p_sim(li[i], sims_i) })
  lag <- as.vector(W %*% z)
  quad <- ifelse(z > 0 & lag > 0, "HH", ifelse(z < 0 & lag < 0, "LL", ifelse(z > 0 & lag < 0, "HL", "LH")))
  list(I = I, p = p_sim(I, sims), local = data.frame(Ii = li, p_sim = lp, quadrant = quad, sig = lp < 0.05))
}
global_rows <- list(); lisa_rows <- list()
for (mx in names(Wfun)) for (yr in CFG$years) {
  W <- Wfun[[mx]](yr); if (is.null(W)) next
  W0 <- if (mx == "hsr") w_hsr_raw(hsr_tau[[hsr_snapshot(yr)]]) else W
  a <- active_W(W0); if (length(a$keep) < 4) next
  x <- dat$D_a[dat$year == yr][a$keep]                      # annual-scale index D_a, ordered by CITY16
  r <- moran_year(x, a$W)
  global_rows[[length(global_rows) + 1]] <- data.frame(matrix = mx, year = yr, N = length(a$keep), I = r$I, p_perm = r$p)
  if (mx %in% c("queen", "hsr"))
    lisa_rows[[length(lisa_rows) + 1]] <- cbind(matrix = mx, year = yr, city_code = CITY16[a$keep], r$local)
}
moran_tab <- do.call(rbind, global_rows); lisa_tab <- do.call(rbind, lisa_rows)
write.csv(moran_tab, outp(tag_s("moran_global.csv")), row.names = FALSE)
write.csv(lisa_tab,  outp(tag_s("lisa_local.csv")),   row.names = FALSE)
# Table A1: Moran's I [permutation p] by matrix, with the thesis values for comparison
A1_THESIS <- data.frame(year = 2011:2023,
  queen = c(-0.274,-0.273,-0.231,-0.257,-0.212,-0.183,-0.204,-0.222,-0.148,-0.058,-0.172,-0.170,-0.171),
  knn4  = c(-0.238,-0.243,-0.228,-0.276,-0.199,-0.174,-0.207,-0.205,-0.146,-0.127,-0.165,-0.183,-0.197),
  idw   = c(-0.107,-0.110,-0.108,-0.109,-0.096,-0.091,-0.100,-0.100,-0.087,-0.074,-0.084,-0.092,-0.079),
  hsr   = c(NA, NA, NA, NA, -0.085,-0.089,-0.113,-0.113,-0.074,-0.055,-0.086,-0.089,-0.081))
cat("\n========== Table A1: Global Moran's I of D_a (3.39)-(3.40), 999 permutations ==========\n")
for (mx in unique(moran_tab$matrix)) {
  r  <- moran_tab[moran_tab$matrix == mx, ]
  th <- A1_THESIS[[mx]][match(r$year, A1_THESIS$year)]
  print(data.frame(matrix = mx, year = r$year, N = r$N, I = round(r$I, 3), p_perm = round(r$p_perm, 3),
                   I_thesis = th, Same_as_thesis = ifelse(abs(round(r$I, 3) - th) < 0.0015, "yes", "NO")),
        row.names = FALSE)
}

# ---- LISA: local Moran's I, queen matrix, Eq (3.41)-(3.42) -----------------
lisa_q <- lisa_tab[lisa_tab$matrix == "queen", ]
lisa_q$city    <- CITY_NAME[as.character(lisa_q$city_code)]
lisa_q$cluster <- ifelse(lisa_q$sig, lisa_q$quadrant, "not significant")

cat("\nLocal Moran's I (LISA), queen matrix, D_a (999 permutations):\n")
for (yr in CFG$years) {
  cat("\n--", yr, "--\n")
  print(lisa_q[lisa_q$year == yr, c("city", "Ii", "p_sim", "quadrant", "cluster")],
        digits = 3, row.names = FALSE)
}
cat("\nSignificant LISA clusters (p < 0.05), number of cities by year:\n")
print(table(lisa_q$year, lisa_q$cluster))

# ======================= SLX exposure terms and models =========================
xw <- function(v) tapply(dat[[v]], list(dat$city_code, dat$year), sum)    # city x year matrix (ordered by CITY16)
exposure <- function(v, Wf, shift = 0) {        # shift = 0 -> (3.43);  shift = 1 -> (3.44): W_{t-1} x_{t-1}
  M <- xw(v); out <- rep(NA_real_, nrow(dat))
  for (yr in unique(dat$year)) {
    s <- yr - shift; if (!as.character(s) %in% colnames(M)) next
    W <- Wf(s); if (is.null(W)) next
    val <- as.vector(W %*% M[, as.character(s)]); val[rowSums(W) == 0] <- NA   # structurally isolated -> NA, not 0
    i <- dat$year == yr; out[i] <- val[match(dat$city_code[i], CITY16)] }
  out
}
slx <- NULL
for (mx in names(Wfun)) {
  d <- dat; d$WD <- exposure("Dg", Wfun[[mx]], 0); d$WD_l1 <- exposure("Dg", Wfun[[mx]], 1)
  for (cv in CONTROLS) d[[paste0("W", cv)]] <- exposure(cv, Wfun[[mx]], 0)
  WX <- paste0("W", CONTROLS)
  specs <- list(BASE = c("Dg", CONTROLS), LITE = c("Dg", "WD", "WD_l1", CONTROLS), FULL = c("Dg", "WD", "WD_l1", CONTROLS, WX))
  if (mx == "hsr") d <- d[d$year >= 2015, ]
  for (yv in c("HQD_ent", "HQD_eq")) for (sp in names(specs)) {     # Table 5.15 / A5: Panel A entropy, Panel B equal weight
    f <- twfe(d, yv, specs[[sp]])
    slx <- rbind(slx, cbind(outcome = yv, matrix = mx, spec = sp,
                            coef_tab(f, paste(yv, mx, sp), intersect(c("Dg", "WD", "WD_l1"), specs[[sp]]))))
  }
}
write.csv(slx, outp(tag_s("results_slx.csv")), row.names = FALSE)
cat("\n========== Table 5.15 / A5: SLX (own-city Dg, WD = exposure, WD_l1 = lagged exposure) ==========\n")
cat("Compare se_k / p_k with the thesis standard errors.\n")
print(slx[, c("outcome","matrix","spec","term","est","se","se_k","p_k","n")], digits = 3, row.names = FALSE)

# ======================= HSR threshold models + multiple testing ================
thr <- NULL
for (yv in c("HQD_ent", "HQD_eq")) for (cc in c(120, 180, 240)) for (side in c("within", "beyond")) {
  Wf <- function(t) { s <- hsr_snapshot(t); if (is.na(s)) return(NULL)
    W0 <- if (side == "within") w_hsr_within(hsr_tau[[s]], cc) else w_hsr_beyond(hsr_tau[[s]], cc); row_standardise(W0) }
  d <- dat; d$WDc <- exposure("Dg", Wf, 0); d <- d[d$year >= 2015, ]
  if (sum(!is.na(d$WDc)) < 30) {                      # too few eligible partners: keep the row, flag it, never drop silently
    warning("Threshold test ", yv, " ", cc, " ", side, " has <30 eligible observations; reported as NA.")
    thr <- rbind(thr, data.frame(outcome = yv, cutoff = cc, side = side, model = paste(yv, cc, side), term = "WDc",
                                 est = NA, se = NA, t = NA, p = NA, se_k = NA, p_k = NA,
                                 n = sum(!is.na(d$WDc)), G = NA, r2_within = NA)); next }
  f <- twfe(d, yv, c("Dg", "WDc", CONTROLS))
  thr <- rbind(thr, cbind(outcome = yv, cutoff = cc, side = side, coef_tab(f, paste(yv, cc, side), "WDc")))
}
thr$p_holm <- NA_real_; thr$p_bh <- NA_real_
for (yv in unique(thr$outcome)) { i <- thr$outcome == yv                 # adjust within each dependent-variable family
  thr$p_holm[i] <- p.adjust(thr$p[i], "holm", n = 6); thr$p_bh[i] <- p.adjust(thr$p[i], "BH", n = 6) }   # family = 3 cutoffs x 2 sides
write.csv(thr, outp(tag_s("results_hsr_thresholds.csv")), row.names = FALSE)
cat("\n========== Tables 5.16-5.17: HSR threshold exposure (3.46)-(3.47), raw / Holm / BH ==========\n")
print(thr[, c("outcome","cutoff","side","est","se","se_k","p","p_holm","p_bh","n")], digits = 3, row.names = FALSE)
