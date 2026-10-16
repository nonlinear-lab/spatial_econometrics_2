# =====================================================================
# R1  Standardisation and Entropy Weighting                (Section 3.4.1)
# ---------------------------------------------------------------------
# Equations implemented
#   (3.1)  min-max standardisation, positive indicator
#   (3.2)  min-max standardisation, negative indicator
#   (3.3)  proportional contribution p_ij
#   (3.4)  entropy value e_j  (k = 1/ln m, 0*ln0 = 0)
#   (3.5)  divergence coefficient g_j = 1 - e_j
#   (3.6)  entropy weight w_j
#   (3.7)  composite score U_i
# Applied to: (a) TPI annual-scale pipeline (m = 16 per year, weights per year)
#             (b) TPI unified-scale pipeline (m = 208, one fixed weight vector)
#             (c) benchmark HQD index + equal-weight, PCA and leave-one-dimension-out
#                 alternatives (needs hqd_controls.xlsx; otherwise DEMO data)
# Outputs: output/tpi_scores.rds|csv, tpi_entropy_weights_*.csv, hqd_indices.rds
# =====================================================================
# Optional: Set your Working Directory (WD) where the graphs will be saved
setwd("C:/R/spatial") 
source("7. R_code for reproducibility.R")
# (removed: CFG$out_dir <- "." -- outputs must go to the "output" folder,
#  because scripts 2-5 read them from there)
options(scipen = 999)

# ---- Equations ------------------------------------------------------
std_positive <- function(a) {                                   # (3.1)
  r <- max(a) - min(a)
  if (r == 0) stop("Indicator has zero range; cannot standardise.")
  (a - min(a)) / r
}
std_negative <- function(a) {                                   # (3.2)
  r <- max(a) - min(a)
  if (r == 0) stop("Indicator has zero range; cannot standardise.")
  (max(a) - a) / r
}
# prop_contrib <- function(X) sweep(X, 2, colSums(X), "/")        # (3.3)  X: m x n matrix
prop_contrib <- function(X) {
  P <- X
  for (j in 1:ncol(X)) { 
    P[, j] <- X[, j] / sum(X[, j])
  }
  return(P)
}
entropy_value <- function(P) {                                  # (3.4)
  m <- nrow(P)
  k <- 1 / log(m)
  e <- numeric(ncol(P))              # empty vector, one slot per indicator
  names(e) <- colnames(P)
  
  for (j in 1:ncol(P)) {             # go through each indicator (column)
    p <- P[, j]
    total <- 0
    for (i in 1:m) {                 # go through each city-year (row)
      if (p[i] > 0) {
        total <- total + p[i] * log(p[i])
      }
    }
    e[j] <- -k * total
  }
  return(e)
}
divergence <- function(e) {
  e1 <- 1 - e                                 # (3.5)
  return(e1)
}
entropy_weight <- function(g) {
  g1 <- g / sum(g)                        # (3.6)
  return(g1)
}
composite_score <- function(X, w) {                              # (3.7)
  n <- nrow(X)
  U <- numeric(n)                    # empty vector, one slot per city-year
  
  for (i in 1:n) {                   # go through each city-year (row)
    total <- 0
    for (j in 1:ncol(X)) {           # go through each indicator (column)
      total <- total + X[i, j] * w[j]
    }
    U[i] <- total
  }
  return(U)
}
# One full pass of (3.1)-(3.7) on a data block A (rows = observations, cols = indicators)
entropy_index <- function(A, direction = rep("+", ncol(A))) {
  A <- as.matrix(A)
  stopifnot(length(direction) == ncol(A))
  X <- sapply(seq_len(ncol(A)), 
      function(j) if (direction[j] == "+") std_positive(A[, j]) 
      else std_negative(A[, j]))
  colnames(X) <- colnames(A)
  P <- prop_contrib(X)
  e <- entropy_value(P)
  g <- divergence(e)
  w <- entropy_weight(g)
  list(X = X, e = e, g = g, w = w, U = composite_score(X, w))
}

# ---- (a)(b) TPI pipelines ---------------------------------------------
panel <- load_tpi_panel()
stopifnot(nrow(panel) == 208)

w_rows <- list(); sc <- panel[, c("city_code", "city", "year")]
for (s in names(TPI_GROUPS)) { sc[[paste0("U", s, "_a")]] <- NA_real_; sc[[paste0("U", s, "_u")]] <- NA_real_ }

# Unified scale: pooled min-max and ONE fixed weight vector (m = 208)
for (s in names(TPI_GROUPS)) {
  ind <- TPI_GROUPS[[s]]; res <- entropy_index(panel[, ind])
  sc[[paste0("U", s, "_u")]] <- res$U
  w_rows[[length(w_rows) + 1]] <- data.frame(pipeline = "unified", year = NA, subsystem = s,
                                             indicator = ind, e = res$e, g = res$g, w = res$w)
}
# Annual scale: year-specific min-max and year-specific weights (m = 16)
for (yr in CFG$years) {
  idx <- which(panel$year == yr); stopifnot(length(idx) == 16)
  for (s in names(TPI_GROUPS)) {
    ind <- TPI_GROUPS[[s]]; res <- entropy_index(panel[idx, ind])
    sc[[paste0("U", s, "_a")]][idx] <- res$U
    w_rows[[length(w_rows) + 1]] <- data.frame(pipeline = "annual", year = yr, subsystem = s,
                                               indicator = ind, e = res$e, g = res$g, w = res$w)
  }
}
wts <- do.call(rbind, w_rows)
write.csv(wts[wts$pipeline == "unified", ], outp("tpi_entropy_weights_unified.csv"), row.names = FALSE)
write.csv(wts[wts$pipeline == "annual", ],  outp("tpi_entropy_weights_annual.csv"),  row.names = FALSE)
saveRDS(sc, outp("tpi_scores.rds")); write.csv(sc, outp("tpi_scores.csv"), row.names = FALSE)
saveRDS(panel, outp("tpi_panel.rds"))
cat("\nUnified-scale entropy weights (3.6):\n"); print(wts[wts$pipeline == "unified", c("subsystem","indicator","e","g","w")], digits = 4, row.names = FALSE)

# ---- (c) HQD index (benchmark entropy, equal-weight, PCA, leave-one-dimension-out) ---
ex <- load_extra(panel)
hq <- merge(panel[, c("city_code", "year")], ex, by = c("city_code", "year"))
hq <- hq[order(hq$city_code, hq$year), ]; stopifnot(nrow(hq) == 208)
hq$Sha1 <- log(hq$Sha1)                         # ln of per-capita public-service expenditure (Sec. 3.2.1)
# The thesis HQD index uses the LOG of the urban-rural income ratio (column "Ln(Coor1)").
# With this line the index reproduces Table 5.1 exactly (mean 0.1874, SD 0.1044);
# with raw Coor1 it does not (mean 0.1850).
if ("Ln(Coor1)" %in% names(hq)) hq$Coor1 <- hq[["Ln(Coor1)"]]

# Benchmark: pooled min-max (m = 208) + entropy weights jointly over the ten indicators
bench <- entropy_index(hq[, HQD_IND$code], HQD_IND$dir)
hq$HQD_ent <- bench$U
write.csv(data.frame(HQD_IND, e = bench$e, g = bench$g, w = bench$w),
          outp(tag("hqd_entropy_weights.csv")), row.names = FALSE)

# Dimension scores: mean of pooled-standardised indicators within each dimension
Xs <- bench$X; dims <- unique(HQD_IND$dim)
DS <- sapply(dims, function(d) rowMeans(Xs[, HQD_IND$code[HQD_IND$dim == d], drop = FALSE]))
colnames(DS) <- dims
hq$HQD_eq <- rowMeans(DS)                       # equal weights of 1/5 on the five dimensions

# PCA alternative: standardise five dimension scores, correlation-matrix PCA, first three PCs,
# combine with relative explained-variance shares, rescale to [0,1]
pc <- prcomp(DS, center = TRUE, scale. = TRUE)
vs <- pc$sdev^2 / sum(pc$sdev^2)
cat("\nPCA: cumulative variance of first 3 PCs =", round(sum(vs[1:3]), 3), "(thesis: approx. 0.80)\n")
comp <- pc$x[, 1:3] %*% (vs[1:3] / sum(vs[1:3]))
hq$HQD_pca <- as.vector((comp - min(comp)) / (max(comp) - min(comp)))

# Leave-one-dimension-out: equal-weight average of the remaining four dimension scores
for (d in dims) hq[[paste0("HQD_lodo_", d)]] <- rowMeans(DS[, setdiff(dims, d), drop = FALSE])
hq <- cbind(hq, setNames(as.data.frame(DS), paste0("DIM_", dims)))
saveRDS(hq, outp(tag("hqd_indices.rds"))); write.csv(hq, outp(tag("hqd_indices.csv")), row.names = FALSE)
cat("\nHQD indices written:", tag("hqd_indices.csv"), "\n")
