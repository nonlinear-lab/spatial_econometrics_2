# =====================================================================
# R2  Coupling Coordination Model                          (Section 3.4.2)
# ---------------------------------------------------------------------
# Equations implemented
#   (3.8)   coupling component  C_it = [ U_T U_P U_I / ((U_T+U_P+U_I)/3)^3 ]^(1/3)
#   (3.9)   development component  T_it = a U_T + b U_P + g U_I
#   (3.10)  a = b = g = 1/3
#   (3.11)  coordination index  D_it = sqrt(C_it * T_it)
# Run on both pipelines:  annual scale -> D_a ;  unified scale -> D_g (written D_u)
# Also: fixed stage classification [0,.3) [.3,.4) [.4,.5) [.5,.6) [.6,.7) [.7,1]
# Outputs: output/tpi_index.rds|csv, tpi_stage_distribution_*.csv, tpi_D_by_year.csv
# =====================================================================
setwd("C:/R/spatial") 
source("7. R_code for reproducibility.R")
need_output("tpi_scores.rds", "R1_standardisation_entropy.R")
sc <- readRDS(outp("tpi_scores.rds"))

# ---- Equations ------------------------------------------------------
coupling_C <- function(UT, UP, UI) {                         # (3.8)
  s <- (UT + UP + UI) / 3
  C <- (UT * UP * UI / s^3)^(1/3)
  C[s == 0] <- 0                                             # all subsystems at their minimum
  return(C)
}
ALPHA <- BETA <- GAMMA <- 1/3                                # (3.10)
development_T <- function(UT, UP, UI) ALPHA * UT + BETA * UP + GAMMA * UI   # (3.9)
coordination_D <- function(C, T) sqrt(C * T)                 # (3.11)

# ---- Apply to both pipelines -------------------------------------------
for (sfx in c("a", "u")) {
  UT <- sc[[paste0("UT_", sfx)]]; UP <- sc[[paste0("UP_", sfx)]]; UI <- sc[[paste0("UI_", sfx)]]
  sc[[paste0("C_", sfx)]] <- coupling_C(UT, UP, UI)
  sc[[paste0("T_", sfx)]] <- development_T(UT, UP, UI)
  sc[[paste0("D_", sfx)]] <- coordination_D(sc[[paste0("C_", sfx)]], sc[[paste0("T_", sfx)]])
}
stopifnot(all(sc$D_a >= 0 & sc$D_a <= 1), all(sc$D_u >= 0 & sc$D_u <= 1))
sc$Dg <- sc$D_u                                              # D_g: unified-scale index used in Chapter 5

# ---- Stage classification (fixed substantive intervals) -----------------
STAGE_BREAKS <- c(0, 0.3, 0.4, 0.5, 0.6, 0.7, Inf)
STAGE_LABELS <- c("Extremely low [0,0.3)", "Low [0.3,0.4)", "Medium-low [0.4,0.5)",
                  "Medium [0.5,0.6)", "Medium-high [0.6,0.7)", "High [0.7,1]")
stage <- function(D) cut(D, STAGE_BREAKS, labels = STAGE_LABELS, right = FALSE)
sc$stage_a <- stage(sc$D_a); sc$stage_u <- stage(sc$D_u)

saveRDS(sc, outp("tpi_index.rds")); write.csv(sc, outp("tpi_index.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(table(sc$year, sc$stage_a)), outp("tpi_stage_distribution_annual.csv"))
write.csv(as.data.frame.matrix(table(sc$year, sc$stage_u)), outp("tpi_stage_distribution_unified.csv"))
by_year <- aggregate(cbind(UT_u, UP_u, UI_u, C_u, T_u, D_a, D_u) ~ year, sc, mean)
write.csv(by_year, outp("tpi_D_by_year.csv"), row.names = FALSE)
cat("\nMean subsystem scores and coordination indices by year:\n"); print(round(by_year, 4), row.names = FALSE)
cat("\nUnified-scale stage distribution (city-years):\n"); print(table(sc$stage_u))
sc <- readRDS(outp("tpi_index.rds"))
vars <- c("D_a", "D_u", "UT_u", "UP_u", "UI_u")
m <- aggregate(sc[vars], list(year = sc$year), mean)
s <- aggregate(sc[vars], list(year = sc$year), sd)
tab <- data.frame(Year = m$year,
                  sapply(vars, function(v) sprintf("%.4f (%.4f)", m[[v]], s[[v]])))
names(tab) <- c("Year", "Annual-scale D", "Unified-scale Dg", "Transport", "Population", "Industry")
print(tab, row.names = FALSE)
write.csv(tab, outp("table_4_1.csv"), row.names = FALSE)
