################################################################################################################
# Habitat suitability model for Solenosmilia variabilis in the Mauritius region
# Seetapah V. and Jogee S.Y. "Predicting deep-sea coral habitat in a data-poor region: a globally calibrated habitat
# model for Solenosmilia variabilis in the Mauritius region". Submitted to Deep-Sea Research Part II.
#
# R code for the analyses, tables and figures reported in the paper. Run it section by section from the repository
# folder. Section 0 must be run first in every session; after that, each section reads what it needs from OUT_DIR,
# so a section can be rerun on its own once the earlier ones have run.
#
# Software. R 4.4.0 with terra 1.7.78, sf 1.0.20, maxnet 0.1.4, ENMeval 2.0.5.2, blockCV 3.2.0, pROC 1.18.5,
# ecospat 4.1.2, predicts 0.1.19, ggplot2 4.0.1 and patchwork. Random seeds are fixed at 42.
#
# Inputs (set INPUT_DIR below). They are not distributed with this file.
#   Seven predictor rasters on one global WGS 84 grid of 0.0051 x 0.0041 degrees (about 560 x 450 m).
#     bathymetry.tif         GEBCO 2023 (m, negative below sea level)
#     slope_degrees.tif      slope derived from GEBCO 2023 with terra::terrain (degrees)
#     aspect_sin.tif         sine of aspect derived from GEBCO 2023
#     aspect_cos.tif         cosine of aspect derived from GEBCO 2023
#     temperature.tif        Bio-ORACLE v3 near-bottom temperature, 2000 to 2020 mean (deg C)
#     kinetic_energy.tif     KE = 0.5 * rho * U^2 from Bio-ORACLE v3 near-bottom current speed U (m/s), rho = 1025 kg/m3 (N/m2)
#     geomorphic_class.tif   Seafloor Geomorphic Features Map (Harris et al. 2014), integer class codes
#   Solenosmilia_variabilis_OBIS_screened.csv   OBIS records after quality screening (8,469 records); columns
#                                               decimalLongitude, decimalLatitude, depth, year
#   Solenosmilia_variabilis_sample400.csv       the 400-record sample used for training (375 with complete predictors);
#                                               same columns
#   OBIS_other_scleractinians.csv               records of the four other scleractinian species (target-group test);
#                                               columns species, lon, lat
#   boundaries/EEZ_MU.shp, Chagos.shp, Tromelin.shp, JMCA.shp
#                                               Mauritian EEZ around Mauritius, Rodrigues, Agalega and St Brandon, around
#                                               the Chagos Archipelago and around Tromelin, and the Joint Management Area
################################################################################################################


# ---- 0. Setup --------------------------------------------------------------------------------------------------
INPUT_DIR <- "inputs"     # folder holding the inputs listed above
OUT_DIR   <- "outputs"    # results are written here
OUT_TAB <- file.path(OUT_DIR, "tables"); OUT_FIG <- file.path(OUT_DIR, "figures"); OUT_R <- file.path(OUT_DIR, "rasters")
for (d in c(OUT_DIR, OUT_TAB, OUT_FIG, OUT_R)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
pkgs <- c("terra", "sf", "maxnet", "ENMeval", "blockCV", "pROC", "ecospat", "predicts", "ggplot2", "patchwork", "viridisLite")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install the missing packages first:\n  install.packages(c(\"", paste(missing, collapse = "\", \""), "\"))")
N_CORES <- min(4, parallel::detectCores())                # cores used by ENMeval; results do not depend on this
suppressPackageStartupMessages({
  library(terra); library(sf); library(maxnet); library(ENMeval); library(blockCV); library(pROC)
  library(ecospat); library(predicts); library(ggplot2); library(patchwork)
})
terraOptions(progress = 0, memfrac = 0.6)
layers <- c(temperature = "temperature.tif", slope_deg = "slope_degrees.tif", kinetic = "kinetic_energy.tif",
            class_cat = "geomorphic_class.tif", bathy = "bathymetry.tif", sin_aspect = "aspect_sin.tif", cos_aspect = "aspect_cos.tif")
vars <- names(layers)                                    # predictor order used by the models
cont <- setdiff(vars, "class_cat")                       # the six continuous predictors
lev  <- as.character(c(1:6, 8:22))                       # geomorphic class codes present in the layer
class_labels <- data.frame(level = lev, label = c("Abyss", "Basins", "Bridges", "Canyons", "Escarpments", "Fans", "Guyots",
  "Hadal", "Plateaus", "Ridges", "Rift valleys", "Rises", "Seamounts", "Shelf", "Shelf valleys", "Sills", "Slope",
  "Spreading ridges", "Terraces", "Trenches", "Troughs"))
# projection window, 52.39 to 75.85 E and 23.81 to 2.29 S, aligned with the predictor grid
window <- ext(52.38948062589422, 52.38948062589422 + 4634 * 0.005063335851830888,
              -2.290475240868327 - 5273 * 0.0040811212187524685, -2.290475240868327)
pred <- function(m, d) as.numeric(predict(m, d, type = "logistic"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))   # versions used for this run


# ---- 1. Presences, background and model table ------------------------------------------------------------------
stk <- rast(file.path(INPUT_DIR, layers)); names(stk) <- vars
vars_mt <- c("temperature", "slope_deg", "kinetic", "bathy", "sin_aspect", "cos_aspect", "class_cat")
as_model_df <- function(v) { v <- v[, vars_mt]; v$class_cat <- factor(as.character(v$class_cat), levels = lev); v }
occ <- read.csv(file.path(INPUT_DIR, "Solenosmilia_variabilis_sample400.csv"))
pv  <- extract(stk, cbind(occ$decimalLongitude, occ$decimalLatitude))
ok  <- complete.cases(pv)
occ <- occ[ok, ]; pv <- as_model_df(pv[ok, ])           # 375 presences
# background, 10,000 random cells with complete predictor values over the global grid
set.seed(42)
bg <- spatSample(stk, 60000, method = "random", na.rm = TRUE, values = TRUE, cells = TRUE, xy = TRUE)
bg <- bg[complete.cases(bg[, vars]), ][1:10000, ]
bv <- as_model_df(bg); keep <- !is.na(bv$class_cat); bv <- bv[keep, ]; bg <- bg[keep, ]
X  <- rbind(pv, bv); Y <- c(rep(1, nrow(pv)), rep(0, nrow(bv)))
xy <- rbind(cbind(occ$decimalLongitude, occ$decimalLatitude), cbind(bg$x, bg$y))
cat("presences", nrow(pv), " background", nrow(bv), "\n")
saveRDS(list(X = X, Y = Y, xy = xy), file.path(OUT_DIR, "model_table_375_10000.rds"))


# ---- 2. Spatial block folds (200 km blocks, five folds) --------------------------------------------------------
tb <- readRDS(file.path(OUT_DIR, "model_table_375_10000.rds")); X <- tb$X; Y <- tb$Y; xy <- tb$xy
pts <- st_as_sf(data.frame(xy, Y = Y), coords = c("X1", "X2"), crs = 4326)
set.seed(42)
sb <- cv_spatial(x = pts, column = "Y", size = 200000, k = 5, selection = "random", iteration = 50,
                 progress = FALSE, report = FALSE, plot = FALSE)
fold <- sb$folds_ids
saveRDS(list(fold = fold, block_m = 200000, seed = 42), file.path(OUT_DIR, "cv_folds_blockCV200km.rds"))
# table in the form used by ENMeval and the later sections
occs <- data.frame(lon = xy[Y == 1, 1], lat = xy[Y == 1, 2], X[Y == 1, vars])
bgd  <- data.frame(lon = xy[Y == 0, 1], lat = xy[Y == 0, 2], X[Y == 0, vars])
saveRDS(list(occs = occs, bg = bgd, occ_grp = fold[Y == 1], bg_grp = fold[Y == 0], levels = lev),
        file.path(OUT_DIR, "model_table_current.rds"))


# ---- 3. Predictor layers cropped to the projection window, and a coarse (x10, about 5 km) stack -------------------
for (nm in vars) {
  out <- file.path(OUT_R, paste0("aoi_", nm, ".tif"))
  if (!file.exists(out)) crop(rast(file.path(INPUT_DIR, layers[[nm]])), window, snap = "near", filename = out, overwrite = TRUE,
                              wopt = list(gdal = "COMPRESS=DEFLATE", datatype = if (nm == "class_cat") "INT2S" else "FLT4S"))
}
cs <- list()
for (nm in vars) cs[[nm]] <- aggregate(rast(file.path(OUT_R, paste0("aoi_", nm, ".tif"))), fact = 10,
                                       fun = if (nm == "class_cat") "modal" else "mean", na.rm = TRUE)
cs <- rast(cs); names(cs) <- vars
writeRaster(cs, file.path(OUT_R, "aoi_predictors_coarse10.tif"), overwrite = TRUE)


# ---- 4. Tuning with ENMeval on the spatial folds and selection by minimum AICc -----------------------------------
set.seed(42)
tb <- readRDS(file.path(OUT_DIR, "model_table_current.rds"))
ev <- ENMevaluate(occs = tb$occs, bg = tb$bg, envs = NULL, algorithm = "maxnet", partitions = "user",
                  user.grp = list(occs.grp = tb$occ_grp, bg.grp = tb$bg_grp),
                  tune.args = list(fc = c("L", "LQ", "LQH", "LQHP", "H"), rm = c(0.5, 1, 1.5, 2, 2.5, 3, 4)),
                  categoricals = "class_cat", doClamp = TRUE, parallel = TRUE, numCores = N_CORES, quiet = TRUE)
saveRDS(ev, file.path(OUT_DIR, "enmeval_object_random_bg.rds"))
r <- ev@results; r <- r[order(r$delta.AICc, -r$auc.val.avg), ]                    # Table S2, selected model first
write.csv(r, file.path(OUT_TAB, "enmeval_tuning_results.csv"), row.names = FALSE)
best <- as.character(r$tune.args[1])                                              # LQH, regularisation multiplier 2
saveRDS(ev@models[[best]], file.path(OUT_DIR, "final_model_current.rds"))
cat("selected:", best, " ncoef", r$ncoef[1], "\n")


# ---- 5. Evaluation, permutation importance, response curves, bivariate surfaces, fold uncertainty (Table 2, Figs 4 to 6, 8)
set.seed(42)
col1 <- viridisLite::viridis(1, begin = 0.35)
tb  <- readRDS(file.path(OUT_DIR, "model_table_current.rds"))
res <- read.csv(file.path(OUT_TAB, "enmeval_tuning_results.csv"))
fc  <- tolower(as.character(res$fc[1])); rm_ <- as.numeric(res$rm[1])
X <- rbind(tb$occs[, vars], tb$bg[, vars]); Y <- c(rep(1, nrow(tb$occs)), rep(0, nrow(tb$bg)))
X$class_cat <- factor(as.character(X$class_cat), levels = lev)
fold <- c(tb$occ_grp, tb$bg_grp)
fm <- maxnet.formula(Y, X, classes = fc)
m_final <- readRDS(file.path(OUT_DIR, "final_model_current.rds"))
# 5.1 spatial-fold evaluation (withheld presences against the background of the same fold)
cv <- data.frame(); fold_models <- list()
for (k in 1:5) {
  tr <- fold != k; te <- fold == k
  m <- maxnet(Y[tr], X[tr, ], f = fm, regmult = rm_); fold_models[[k]] <- m
  p <- pred(m, X[te, ]); rr <- roc(Y[te], p, quiet = TRUE)
  cv <- rbind(cv, data.frame(fold = k, n_test_presences = sum(Y[te] == 1), n_test_background = sum(Y[te] == 0),
                             auc = as.numeric(auc(rr)), max_tss = max(rr$sensitivities + rr$specificities - 1),
                             boyce = ecospat.boyce(fit = p[Y[te] == 0], obs = p[Y[te] == 1], nclass = 0, PEplot = FALSE)$cor))
}
p_all <- pred(m_final, X); auc_train <- as.numeric(auc(roc(Y, p_all, quiet = TRUE)))
summ <- data.frame(metric = c("AUC", "TSS (max)", "Continuous Boyce index", "AUC full model on training data"),
                   mean = c(mean(cv$auc), mean(cv$max_tss), mean(cv$boyce), auc_train), sd = c(sd(cv$auc), sd(cv$max_tss), sd(cv$boyce), NA))
write.csv(cv, file.path(OUT_TAB, "v2_cv_per_fold.csv"), row.names = FALSE)
write.csv(summ, file.path(OUT_TAB, "v2_cv_summary.csv"), row.names = FALSE)
saveRDS(fold_models, file.path(OUT_DIR, "v2_fold_models.rds"))
# 5.2 permutation importance (decrease in AUC, 10 permutations)
imp <- data.frame()
for (v in vars) {
  d <- replicate(10, { Xp <- X; Xp[[v]] <- sample(Xp[[v]]); auc_train - as.numeric(auc(roc(Y, pred(m_final, Xp), quiet = TRUE))) })
  imp <- rbind(imp, data.frame(variable = v, auc_full = auc_train, auc_drop_mean = mean(d), auc_drop_sd = sd(d)))
}
imp$relative_importance_pct <- 100 * pmax(imp$auc_drop_mean, 0) / sum(pmax(imp$auc_drop_mean, 0))
imp <- imp[order(-imp$auc_drop_mean), ]
write.csv(imp, file.path(OUT_TAB, "v2_permutation_importance.csv"), row.names = FALSE)
# 5.3 response curves, one predictor over its 0.1 to 99.9 percentile range, the others at their mean over the training
#     records (modal class for the geomorphic layer); categorical response; effect sizes
P <- X[Y == 1, ]
base <- as.data.frame(lapply(vars, function(v) if (v == "class_cat") factor(names(which.max(table(P$class_cat))), levels = lev) else mean(P[[v]])))
names(base) <- vars
pdp1 <- function(v, xs) { d <- base[rep(1, length(xs)), ]; d[[v]] <- if (v == "class_cat") factor(xs, levels = lev) else xs; pred(m_final, d) }
rc <- list(); eff <- data.frame()
for (v in cont) {
  xs <- seq(quantile(X[[v]], 0.001), quantile(X[[v]], 0.999), length.out = 200)
  y <- pdp1(v, xs); rc[[v]] <- data.frame(x = xs, y = y)
  write.csv(rc[[v]], file.path(OUT_TAB, paste0("v2_response_", v, ".csv")), row.names = FALSE)
  q <- xs[y >= 0.9 * max(y)]
  eff <- rbind(eff, data.frame(variable = v, suit_min = min(y), suit_max = max(y), effect_size = max(y) - min(y),
                               x_at_max = xs[which.max(y)], x90_low = min(q), x90_high = max(q)))
}
cat_resp <- data.frame(level = lev, suitability = pdp1("class_cat", lev), n_presences = as.integer(table(factor(P$class_cat, levels = lev))))
cat_resp <- merge(cat_resp, class_labels, by = "level", all.x = TRUE)
eff <- rbind(eff, data.frame(variable = "class_cat", suit_min = min(cat_resp$suitability), suit_max = max(cat_resp$suitability),
                             effect_size = diff(range(cat_resp$suitability)), x_at_max = NA, x90_low = NA, x90_high = NA))
write.csv(cat_resp[order(-cat_resp$suitability), ], file.path(OUT_TAB, "v2_response_class_cat.csv"), row.names = FALSE)
write.csv(eff, file.path(OUT_TAB, "v2_effect_sizes.csv"), row.names = FALSE)
# 5.4 Figures 4 and 5
lab <- c(bathy = "Bathymetry", temperature = "Temperature", class_cat = "Geomorphic feature", slope_deg = "Seafloor slope",
         cos_aspect = "Aspect (cosine)", sin_aspect = "Aspect (sine)", kinetic = "Kinetic energy")
imp$label <- factor(lab[imp$variable], levels = rev(lab[imp$variable]))
p4 <- ggplot(imp, aes(x = auc_drop_mean, y = label)) + geom_col(fill = col1) +
  geom_errorbar(aes(xmin = auc_drop_mean - auc_drop_sd, xmax = auc_drop_mean + auc_drop_sd), width = 0.25) +
  geom_text(aes(x = auc_drop_mean + auc_drop_sd, label = sprintf("%.1f%%", relative_importance_pct)), hjust = -0.3, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Decrease in AUC when the predictor is permuted (mean and SD of 10 permutations)", y = NULL) + theme_minimal(base_size = 12)
ggsave(file.path(OUT_FIG, "Fig4_permutation_importance.png"), p4, width = 7.5, height = 4, dpi = 300)
ymax <- max(sapply(rc, function(d) max(d$y)), cat_resp$suitability) * 1.05
panel <- function(v, xlab, tag) ggplot(rc[[v]], aes(x, y)) + geom_line(colour = col1, linewidth = 0.9) +
  scale_y_continuous(limits = c(0, ymax)) + labs(x = xlab, y = "Suitability", tag = tag) + theme_minimal(base_size = 11)
cat_resp$label <- factor(cat_resp$label, levels = cat_resp$label[order(cat_resp$suitability)])
pg <- ggplot(cat_resp, aes(x = suitability, y = label)) + geom_col(fill = col1) +
  geom_text(aes(label = ifelse(n_presences == 0, "", paste0("n=", n_presences))), hjust = -0.1, size = 2.8) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) + labs(x = "Suitability", y = NULL, tag = "g") + theme_minimal(base_size = 11)
fig5 <- (panel("bathy", "Bathymetry (m)", "a") | panel("temperature", "Temperature (\u00b0C)", "b") | panel("slope_deg", "Seafloor slope (degrees)", "c")) /
        (panel("kinetic", expression(Kinetic~energy~(N~m^-2)), "d") | panel("cos_aspect", "Aspect (cosine)", "e") | panel("sin_aspect", "Aspect (sine)", "f")) /
        pg + plot_layout(heights = c(1, 1, 1.6))
ggsave(file.path(OUT_FIG, "Fig5_response_curves.png"), fig5, width = 10, height = 11, dpi = 300)
# 5.5 bivariate partial dependence (Figure 6)
labu <- c(bathy = "Bathymetry (m)", temperature = "Temperature (\u00b0C)", slope_deg = "Seafloor slope (degrees)", kinetic = expression(Kinetic~energy~(N~m^-2)))
pairs <- list(c("bathy", "slope_deg"), c("bathy", "temperature"), c("bathy", "kinetic"), c("slope_deg", "kinetic"))
pdp <- list(); tags <- c("a", "b", "c", "d")
for (i in seq_along(pairs)) {
  v1 <- pairs[[i]][1]; v2 <- pairs[[i]][2]
  g <- expand.grid(x1 = seq(quantile(X[[v1]], .005), quantile(X[[v1]], .995), length.out = 60),
                   x2 = seq(quantile(X[[v2]], .005), quantile(X[[v2]], .995), length.out = 60))
  d <- base[rep(1, nrow(g)), ]; d[[v1]] <- g$x1; d[[v2]] <- g$x2; g$y <- pred(m_final, d)
  write.csv(g, file.path(OUT_TAB, sprintf("v2_pdp2d_%s_x_%s.csv", v1, v2)), row.names = FALSE)
  pdp[[i]] <- ggplot(g, aes(x1, x2, fill = y)) + geom_raster() + scale_fill_viridis_c(name = "Suitability", limits = c(0, ymax)) +
    labs(x = labu[v1], y = labu[v2], tag = tags[i]) + theme_minimal(base_size = 10)
}
fig6 <- (pdp[[1]] | pdp[[2]]) / (pdp[[3]] | pdp[[4]]) + plot_layout(guides = "collect")
ggsave(file.path(OUT_FIG, "Fig6_bivariate_pdp.png"), fig6, width = 9, height = 7.5, dpi = 300)
# 5.6 standard deviation across the five fold models on the coarse grid (Figure 8)
cs <- rast(file.path(OUT_R, "aoi_predictors_coarse10.tif"))
cells <- which(!is.na(values(cs[[1]]))); D <- as.data.frame(cs[cells])[, vars]
ok <- complete.cases(D) & as.character(D$class_cat) %in% lev; D$class_cat <- factor(as.character(D$class_cat), levels = lev)
Pm <- sapply(fold_models, function(m) pred(m, D[ok, ]))
mk <- function(v) { r <- rast(cs[[1]]); values(r) <- NA; r[cells[ok]] <- v; r }
sd_r <- mk(apply(Pm, 1, sd)); names(sd_r) <- "sd_folds"
mean_r <- mk(rowMeans(Pm)); names(mean_r) <- "mean_folds"
full_c <- mk(pred(m_final, D[ok, ])); names(full_c) <- "final_model"
writeRaster(c(mean_r, sd_r, full_c), file.path(OUT_R, "v2_coarse10_fold_mean_sd_final.tif"), overwrite = TRUE)
to_df <- function(r) { d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"; d }
p8 <- ggplot(to_df(sd_r), aes(x, y, fill = v)) + geom_raster() + coord_sf(crs = 4326, expand = FALSE) +
  scale_fill_viridis_c(option = "cividis", name = "SD") + labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) + theme(panel.background = element_rect(fill = "grey95", colour = NA))
ggsave(file.path(OUT_FIG, "Fig8_uncertainty_sd_coarse10.png"), p8, width = 8, height = 7, dpi = 300)


# ---- 6. Extrapolation diagnostics, MESS and ExDet (Figure 9) ---------------------------------------------------
cs <- rast(file.path(OUT_R, "aoi_predictors_coarse10.tif"))
tb <- readRDS(file.path(OUT_DIR, "model_table_current.rds"))
refP <- tb$occs[, cont]; refPB <- rbind(tb$occs[, cont], tb$bg[, cont])
m1 <- predicts::mess(cs[[cont]], refP, full = FALSE);  names(m1) <- "mess_presences"
m2 <- predicts::mess(cs[[cont]], refPB, full = FALSE); names(m2) <- "mess_presences_background"
cells <- which(!is.na(values(cs[[1]]))); D <- as.data.frame(cs[cells])[, cont]; ok <- complete.cases(D)
set.seed(42); refS <- refPB[sample(nrow(refPB), 4000), ]          # subsample of the reference for the Mahalanobis step
cmn <- ecospat.climan(refS, D[ok, ])                              # ExDet, NT1 < 0 univariate and NT2 > 1 combinational novelty
ex <- rast(cs[[1]]); values(ex) <- NA; ex[cells[ok]] <- cmn; names(ex) <- "exdet"
writeRaster(c(m1, m2, ex), file.path(OUT_R, "v2_coarse10_mess_exdet.tif"), overwrite = TRUE)
fin <- rast(file.path(OUT_R, "v2_coarse10_fold_mean_sd_final.tif"))[["final_model"]]
cls <- classify(fin, cbind(c(-Inf, 0.3, 0.6), c(0.3, 0.6, Inf), 1:3), right = FALSE)
ca <- cellSize(cls, unit = "km", mask = TRUE)
share <- function(flag) { z <- zonal(ca * flag, cls, "sum", na.rm = TRUE)[, 2]; a <- zonal(ca, cls, "sum", na.rm = TRUE)[, 2]; round(100 * z / a, 2) }
tab <- data.frame(class = c("low (<0.3)", "medium (0.3-0.6)", "high (>0.6)"),
                  mess_presences_negative_pct = share(m1 < 0), mess_presences_background_negative_pct = share(m2 < 0),
                  exdet_type1_pct = share(ex < 0), exdet_type2_pct = share(ex > 1))
write.csv(tab, file.path(OUT_TAB, "v2_extrapolation_by_class_coarse10.csv"), row.names = FALSE)
mdv <- predicts::mess(cs[[cont]], refP, full = TRUE)              # most dissimilar variable relative to the presences
md  <- which.min(mdv[[seq_along(cont)]]); neg <- m1 < 0
f <- as.data.frame(freq(mask(md, neg, maskvalues = FALSE))); f$variable <- cont[f$value]; f$pct_of_negative <- 100 * f$count / sum(f$count)
write.csv(f[, c("variable", "pct_of_negative")], file.path(OUT_TAB, "mess_most_dissimilar_variable.csv"), row.names = FALSE)
to_df2 <- function(r, lo = NULL) { d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"; if (!is.null(lo)) d$v <- pmax(d$v, lo); d }
mapb <- function(d, lab, tag, pal) ggplot(d, aes(x, y, fill = v)) + geom_raster() + coord_sf(crs = 4326, expand = FALSE) + pal +
  labs(x = "Longitude", y = "Latitude", fill = lab, tag = tag) + theme_minimal(base_size = 10) + theme(panel.background = element_rect(fill = "grey95", colour = NA))
d1 <- to_df2(m1, -100); d2 <- to_df2(m2, -100); d3 <- to_df2(ex); d3$v <- pmin(d3$v, 3)
p9a <- mapb(d1, "MESS", "a", scale_fill_viridis_c(limits = c(-100, max(d1$v, d2$v))))
p9b <- mapb(d2, "MESS", "b", scale_fill_viridis_c(limits = c(-100, max(d1$v, d2$v))))
p9c <- mapb(d3, "ExDet", "c", scale_fill_viridis_c(option = "cividis"))
ggsave(file.path(OUT_FIG, "Fig9_mess_exdet_coarse10.png"), (p9a | p9b) / (p9c | plot_spacer()), width = 11, height = 10, dpi = 300)


# ---- 7. Projection at full resolution and area of each suitability class (Table 3) ------------------------------
EEZ_DIR <- file.path(INPUT_DIR, "boundaries")                     # EEZ_MU, Chagos, Tromelin and JMCA shapefiles
tb  <- readRDS(file.path(OUT_DIR, "model_table_current.rds"))
stk <- rast(file.path(OUT_R, paste0("aoi_", vars, ".tif"))); names(stk) <- vars
m_final <- readRDS(file.path(OUT_DIR, "final_model_current.rds"))
wrap <- function(model, data, ...) {
  data$class_cat <- factor(as.character(data$class_cat), levels = lev)
  out <- rep(NA_real_, nrow(data)); ok <- complete.cases(data)
  if (any(ok)) out[ok] <- as.numeric(predict(model, data[ok, ], type = "logistic")); out
}
full <- file.path(OUT_R, "v2_pred_full_resolution.tif")
pr <- predict(stk, m_final, fun = wrap, na.rm = FALSE, filename = full, overwrite = TRUE, wopt = list(gdal = "COMPRESS=DEFLATE"))
names(pr) <- "suitability"
shp <- function(n) { v <- vect(file.path(EEZ_DIR, paste0(n, ".shp"))); crs(v) <- "EPSG:4326"; values(v) <- NULL; v }
zones <- list(Chagos = shp("Chagos"), EEZ_MU = shp("EEZ_MU"), Tromelin = shp("Tromelin"), JMA = shp("JMCA"))
mess_p <- rast(file.path(OUT_R, "v2_coarse10_mess_exdet.tif"))[["mess_presences"]]
class_area <- function(v, name) {
  r <- mask(crop(pr, v), v)
  ca <- cellSize(r, unit = "km", mask = TRUE)
  cl <- classify(r, cbind(c(-Inf, 0.3, 0.6), c(0.3, 0.6, Inf), 1:3), right = FALSE)
  neg <- resample(mess_p, r, method = "near") < 0
  a <- zonal(ca, cl, "sum", na.rm = TRUE); n <- zonal(ca * neg, cl, "sum", na.rm = TRUE)
  data.frame(area = name, class = c("low (<0.3)", "medium (0.3-0.6)", "high (>0.6)")[a[, 1]], area_km2 = a[, 2],
             pct_of_area = 100 * a[, 2] / sum(a[, 2]), mess_negative_pct_of_class = 100 * n[, 2] / a[, 2])
}
eez_all <- aggregate(rbind(zones$Chagos, zones$EEZ_MU, zones$Tromelin))
out <- rbind(class_area(zones$Chagos, "Chagos Archipelago"),
             class_area(aggregate(rbind(zones$EEZ_MU, zones$Tromelin)), "Mauritius, Rodrigues, Agal\u00e9ga, St Brandon and Tromelin"),
             class_area(eez_all, "Mauritian EEZ, total"),
             class_area(zones$JMA, "Joint Management Area"))
write.csv(out, file.path(OUT_TAB, "v3_suitability_class_areas_by_subarea.csv"), row.names = FALSE)


# ---- 8. Sensitivity analyses (Tables S3 and S4) ----------------------------------------------------------------
# 8.1 thinning to one record per 5 km, and prediction of the Indian Ocean records by a model fitted without them
set.seed(42)
tb  <- readRDS(file.path(OUT_DIR, "model_table_current.rds"))
res <- read.csv(file.path(OUT_TAB, "enmeval_tuning_results.csv")); fc <- tolower(as.character(res$fc[1])); rm_ <- as.numeric(res$rm[1])
X <- rbind(tb$occs[, vars], tb$bg[, vars]); X$class_cat <- factor(as.character(X$class_cat), levels = lev)
Y <- c(rep(1, nrow(tb$occs)), rep(0, nrow(tb$bg))); fold <- c(tb$occ_grp, tb$bg_grp)
lon <- c(tb$occs$lon, tb$bg$lon); lat <- c(tb$occs$lat, tb$bg$lat)
evalk <- function(p, y) { rr <- roc(y, p, quiet = TRUE); c(auc = as.numeric(auc(rr)), tss = max(rr$sensitivities + rr$specificities - 1),
  boyce = ecospat.boyce(fit = p[y == 0], obs = p[y == 1], nclass = 0, PEplot = FALSE)$cor) }
m_final <- readRDS(file.path(OUT_DIR, "final_model_current.rds"))
fm <- maxnet.formula(Y, X, classes = fc)
np <- nrow(tb$occs); v <- vect(cbind(tb$occs$lon, tb$occs$lat), crs = "EPSG:4326"); dm <- as.matrix(distance(v)) / 1000
keep <- rep(TRUE, np); for (i in seq_len(np)) if (keep[i]) keep[which(dm[i, ] < 5 & seq_len(np) > i)] <- FALSE
sel <- c(which(keep), np + seq_len(nrow(tb$bg)))
Xt <- X[sel, ]; Yt <- Y[sel]; ft <- fold[sel]
cvt <- t(sapply(1:5, function(k) { m <- maxnet(Yt[ft != k], Xt[ft != k, ], f = maxnet.formula(Yt, Xt, classes = fc), regmult = rm_); evalk(pred(m, Xt[ft == k, ]), Yt[ft == k]) }))
mt <- maxnet(Yt, Xt, f = maxnet.formula(Yt, Xt, classes = fc), regmult = rm_)
P <- Xt[Yt == 1, ]; base <- as.data.frame(lapply(vars, function(z) if (z == "class_cat") factor(names(which.max(table(P$class_cat))), levels = lev) else mean(P[[z]]))); names(base) <- vars
curve <- function(m, z, xs) { d <- base[rep(1, length(xs)), ]; d[[z]] <- xs; pred(m, d) }
xs <- seq(-5000, 0, by = 25); yb <- curve(mt, "bathy", xs); xt <- seq(0, 20, by = 0.1); yt <- curve(mt, "temperature", xt)
cs <- rast(file.path(OUT_R, "aoi_predictors_coarse10.tif")); D <- as.data.frame(spatSample(cs, 20000, method = "random", na.rm = TRUE))
D <- D[D$class_cat %in% as.numeric(lev), vars]; D$class_cat <- factor(as.character(D$class_cat), levels = lev)
thin <- data.frame(test = "thinned to one record per 5 km", n_presences = sum(keep),
                   auc = mean(cvt[, "auc"]), auc_sd = sd(cvt[, "auc"]), tss = mean(cvt[, "tss"]), boyce = mean(cvt[, "boyce"]), boyce_sd = sd(cvt[, "boyce"]),
                   depth_peak_m = xs[which.max(yb)], temp_peak_c = xt[which.max(yt)], r_with_final_model = cor(pred(mt, D), pred(m_final, D)))
io <- lon > 20 & lon < 147 & lat < 30                             # Indian Ocean records and background
m_io <- maxnet(Y[!io], X[!io, ], f = fm, regmult = rm_)
e_io <- evalk(pred(m_io, X[io, ]), Y[io])
lobo <- data.frame(test = "Indian Ocean withheld (trained on Atlantic and Pacific)", n_presences = sum(io & Y == 1),
                   auc = e_io["auc"], auc_sd = NA, tss = e_io["tss"], boyce = e_io["boyce"], boyce_sd = NA,
                   depth_peak_m = NA, temp_peak_c = NA, r_with_final_model = cor(pred(m_io, D), pred(m_final, D)))
write.csv(rbind(thin, lobo), file.path(OUT_TAB, "v2_sensitivity_thinning_basin.csv"), row.names = FALSE)

# 8.2 target-group background (Table S3); records of the four other scleractinian species, 10,000 points drawn
#     across their grid cells in proportion to record density, same presences, folds rebuilt, same tuning grid
set.seed(42)
stk <- rast(file.path(INPUT_DIR, layers)); names(stk) <- vars
t0  <- readRDS(file.path(OUT_DIR, "model_table_375_10000.rds"))
P   <- t0$X[t0$Y == 1, ]; Pxy <- t0$xy[t0$Y == 1, ]
tg  <- read.csv(file.path(INPUT_DIR, "OBIS_other_scleractinians.csv"))     # columns species, lon, lat
tg  <- tg[!grepl("Solenosmilia", tg$species) & is.finite(tg$lon) & is.finite(tg$lat), ]
cellsg <- cellFromXY(stk[[1]], cbind(tg$lon, tg$lat))
cnt <- table(cellsg); cell_ids <- as.numeric(names(cnt)); w <- as.numeric(cnt)
uxy <- xyFromCell(stk[[1]], cell_ids); uv <- extract(stk, uxy)[, vars]; okg <- complete.cases(uv)
cell_ids <- cell_ids[okg]; w <- w[okg]; uxy <- uxy[okg, ]; uv <- uv[okg, ]
draw <- sample(seq_along(cell_ids), 10000, replace = TRUE, prob = w)
bxy <- uxy[draw, ]; bvg <- uv[draw, ]
write.csv(data.frame(records = nrow(tg), species = length(unique(tg$species)), cells_with_complete_predictors = length(cell_ids),
                     background_points = nrow(bvg), unique_cells_in_background = length(unique(draw))),
          file.path(OUT_TAB, "target_group_background_summary.csv"), row.names = FALSE)
levg <- sort(unique(c(P$class_cat, bvg$class_cat)))
P$class_cat <- factor(as.character(P$class_cat), levels = as.character(levg)); bvg$class_cat <- factor(as.character(bvg$class_cat), levels = as.character(levg))
xy_all <- rbind(unname(as.matrix(Pxy)), unname(as.matrix(bxy))); colnames(xy_all) <- c("lon", "lat")
ptsg <- st_as_sf(data.frame(xy_all, Y = c(rep(1, nrow(P)), rep(0, nrow(bvg)))), coords = c("lon", "lat"), crs = 4326)
sbg <- cv_spatial(x = ptsg, column = "Y", size = 200000, k = 5, selection = "random", iteration = 50, progress = FALSE, report = FALSE, plot = FALSE)
grp <- sbg$folds_ids
occs_g <- data.frame(lon = Pxy[, 1], lat = Pxy[, 2], P); bg_g <- data.frame(lon = bxy[, 1], lat = bxy[, 2], bvg)
ev_g <- ENMevaluate(occs = occs_g, bg = bg_g, envs = NULL, algorithm = "maxnet", partitions = "user",
                    user.grp = list(occs.grp = grp[seq_len(nrow(P))], bg.grp = grp[-seq_len(nrow(P))]),
                    tune.args = list(fc = c("L", "LQ", "LQH", "LQHP", "H"), rm = c(0.5, 1, 1.5, 2, 2.5, 3, 4)),
                    categoricals = "class_cat", doClamp = TRUE, parallel = TRUE, numCores = N_CORES, quiet = TRUE)
rg <- ev_g@results; rg <- rg[order(rg$or.10p.avg, -rg$auc.val.avg), ]              # sequential criterion (Kass et al. 2021)
write.csv(rg, file.path(OUT_TAB, "enmeval_tuning_results_target_group.csv"), row.names = FALSE)
saveRDS(ev_g@models[[as.character(rg$tune.args[1])]], file.path(OUT_DIR, "final_model_target_group.rds"))
# The selected target-group model (LQHP, regularisation multiplier 4) was evaluated and projected with the code of Sections 5 and 7.


# ---- 9. Predictor correlations (Table S1) and training-presence thresholds for the suitability classes (Section 2.6)
t0 <- readRDS(file.path(OUT_DIR, "model_table_375_10000.rds"))
write.csv(round(cor(t0$X[, cont]), 3), file.path(OUT_TAB, "predictor_correlations.csv"))
tb <- readRDS(file.path(OUT_DIR, "model_table_current.rds")); m_final <- readRDS(file.path(OUT_DIR, "final_model_current.rds"))
Pp <- tb$occs[, vars]; Pp$class_cat <- factor(as.character(Pp$class_cat), levels = lev)
q <- quantile(pred(m_final, Pp), c(0, 0.05, 0.10, 0.5))                  # 5% = 0.27 and 10% = 0.36 in the paper
write.csv(data.frame(quantile = names(q), value = as.numeric(q)), file.path(OUT_TAB, "v2_training_presence_thresholds.csv"), row.names = FALSE)


# ---- 10. Occurrence records, counts, depth and years (Section 2.2, Figures S1 and S2) and nearest-neighbour distances
clean <- read.csv(file.path(INPUT_DIR, "Solenosmilia_variabilis_OBIS_screened.csv"))
s400  <- read.csv(file.path(INPUT_DIR, "Solenosmilia_variabilis_sample400.csv"))
stk   <- rast(file.path(INPUT_DIR, layers)); names(stk) <- vars
s375  <- s400[complete.cases(extract(stk, cbind(s400$decimalLongitude, s400$decimalLatitude))), ]
summ_occ <- function(d, label) {
  dep <- suppressWarnings(as.numeric(d$depth)); yr <- suppressWarnings(as.numeric(d$year))
  data.frame(dataset = label, n_records = nrow(d), depth_min_m = min(dep, na.rm = TRUE), depth_q50_m = median(dep, na.rm = TRUE),
             depth_max_m = max(dep, na.rm = TRUE), n_with_year = sum(!is.na(yr)),
             year_min = min(yr, na.rm = TRUE), year_max = max(yr, na.rm = TRUE))
}
write.csv(rbind(summ_occ(clean, "screened OBIS table"), summ_occ(s400, "sample of 400"), summ_occ(s375, "used by the model")),
          file.path(OUT_TAB, "occurrence_summary.csv"), row.names = FALSE)
s375$depth <- as.numeric(s375$depth)
ggsave(file.path(OUT_FIG, "FigS1_depth_histogram.png"),
       ggplot(s375[!is.na(s375$depth), ], aes(depth)) + geom_histogram(binwidth = 100, fill = viridisLite::viridis(1, begin = 0.35), colour = "white") +
         labs(x = "Depth (m)", y = "Number of records") + theme_minimal(base_size = 12), width = 7, height = 4.5, dpi = 300)
yr <- suppressWarnings(as.numeric(s375$year))
ggsave(file.path(OUT_FIG, "FigS2_year_histogram.png"),
       ggplot(data.frame(year = yr[!is.na(yr)]), aes(year)) + geom_histogram(binwidth = 5, fill = viridisLite::viridis(1, begin = 0.35), colour = "white") +
         labs(x = "Year of record", y = "Number of records") + theme_minimal(base_size = 12), width = 7, height = 4.5, dpi = 300)
nn_km <- function(d) { v <- vect(cbind(d$decimalLongitude, d$decimalLatitude), crs = "EPSG:4326"); dm <- as.matrix(distance(v)) / 1000; diag(dm) <- Inf; apply(dm, 1, min) }
nn <- nn_km(s400)
write.csv(data.frame(set = "sample of 400", nn_median_km = median(nn), n_within_5km = sum(nn < 5)),
          file.path(OUT_TAB, "thinning_check_nearest_neighbour.csv"), row.names = FALSE)


# ---- 11. Figure 3, predictor layers over the study area (colour range = 1st to 99th percentile within the area) ----
cs <- rast(file.path(OUT_R, "aoi_predictors_coarse10.tif"))
to_df <- function(r) { d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"; d }
pan <- function(nm, lab3, tag, q = c(0.01, 0.99)) {
  d <- to_df(cs[[nm]]); lim <- quantile(d$v, q); d$v <- pmin(pmax(d$v, lim[1]), lim[2])
  ggplot(d, aes(x, y, fill = v)) + geom_raster() + scale_fill_viridis_c(name = lab3) + coord_sf(crs = 4326, expand = FALSE) +
    labs(x = NULL, y = NULL, tag = tag) + theme_minimal(base_size = 9) + theme(legend.position = "right", legend.key.height = unit(0.5, "cm"))
}
bd <- to_df(cs[["bathy"]]); bd$v <- pmin(bd$v, 0)
pa <- ggplot(bd, aes(x, y, fill = v)) + geom_raster() + scale_fill_viridis_c(name = "Bathymetry (m)", option = "cividis") +
  coord_sf(crs = 4326, expand = FALSE) + labs(x = NULL, y = NULL, tag = "a") + theme_minimal(base_size = 9) + theme(legend.key.height = unit(0.5, "cm"))
f3 <- (pa | pan("slope_deg", "Slope (degrees)", "b")) / (pan("temperature", "Temperature (\u00b0C)", "c") | pan("kinetic", expression(atop("Kinetic energy", (N~m^-2))), "d")) /
      (pan("cos_aspect", "Aspect, cosine", "e", c(0, 1)) | pan("sin_aspect", "Aspect, sine", "f", c(0, 1)))
ggsave(file.path(OUT_FIG, "Fig3_predictors.png"), f3, width = 9, height = 12, dpi = 300)


# ---- 12. Figures 1 and 7, maps with the maritime boundaries and place names ------------------------------------
EEZ_DIR <- file.path(INPUT_DIR, "boundaries")
zone_names <- c(EEZ_MU = "EEZ, Mauritius, Rodrigues, Agal\u00e9ga and St Brandon", Chagos = "EEZ, Chagos Archipelago",
                Tromelin = "EEZ, Tromelin", JMCA = "Joint Management Area (Mauritius and Seychelles)")
zones <- do.call(rbind, lapply(names(zone_names), function(n) {
  v <- st_read(file.path(EEZ_DIR, paste0(n, ".shp")), quiet = TRUE); st_crs(v) <- 4326
  st_sf(name = n, geometry = st_geometry(st_union(v)))
}))
zones$label <- factor(zone_names[zones$name], levels = zone_names)
okabe <- setNames(c("#000000", "#0072B2", "#CC79A7", "#D55E00"), zone_names)       # Okabe-Ito
ltys  <- setNames(c("solid", "solid", "dotted", "dashed"), zone_names)
zone_title <- "Maritime zones of the Republic of Mauritius"
# place names (decimal degrees, WGS 84); dx and dy move each label away from its marker
isl <- data.frame(
  name = c("Mauritius", "R\u00e9union", "Rodrigues", "Agal\u00e9ga", "St Brandon", "Tromelin", "Seychelles (Mah\u00e9)",
           "Diego Garcia", "Peros Banhos", "Salomon Islands", "Egmont Islands"),
  x = c(57.55, 55.53, 63.42, 56.65, 59.62, 54.52, 55.45, 72.42, 71.82, 72.25, 71.35),
  y = c(-20.25, -21.12, -19.72, -10.42, -16.45, -15.89, -4.62, -7.32, -5.30, -5.33, -6.66))
banks <- data.frame(name = c("Saya de Malha Bank", "Nazareth Bank", "Great Chagos Bank"), x = c(61.3, 60.5, 72.0), y = c(-10.4, -14.6, -6.2))
nudge <- function(d, dx, dy) { d$dx <- dx; d$dy <- dy; d }
lab_box <- function(d, size, fontface = "plain", colour = "black")
  geom_label(data = d, aes(x = x + dx, y = y + dy, label = name), size = size, fontface = fontface, colour = colour,
             fill = scales::alpha("white", 0.88), label.size = 0, label.padding = unit(0.1, "lines"), lineheight = 0.9)
marker <- function(d) geom_point(data = d, aes(x, y), shape = 21, fill = "white", colour = "black", size = 1.6, stroke = 0.5)
leader <- function(d, f = 0.8) geom_segment(data = d, aes(x = x, y = y, xend = x + dx * f, yend = y + dy * f), linewidth = 0.25, colour = "grey15")
lonlab <- function(x) paste0(abs(x), "\u00b0", ifelse(x < 0, "W", "E"))
latlab <- function(y) paste0(abs(y), "\u00b0", ifelse(y < 0, "S", "N"))
deg_axes <- list(scale_x_continuous(labels = lonlab), scale_y_continuous(labels = latlab))
dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
# Figure 1, study area over GEBCO bathymetry
cs <- rast(file.path(OUT_R, "aoi_predictors_coarse10.tif"))
b <- as.data.frame(cs[["bathy"]], xy = TRUE, na.rm = TRUE); names(b)[3] <- "v"; b$v <- pmin(b$v, 0)
i1 <- nudge(isl[isl$name %in% c("Mauritius", "R\u00e9union", "Rodrigues", "Agal\u00e9ga", "St Brandon", "Tromelin",
                                "Seychelles (Mah\u00e9)", "Diego Garcia"), ],
            dx = c(1.9, -1.6, 0, 1.6, -1.6, 0, 1.9, 0), dy = c(0.4, -0.6, 0.85, 0.3, 0.3, 0.85, 0.5, -0.8))
b1 <- nudge(banks, 0, 0)
jma1 <- data.frame(name = "Joint Management\nArea", x = 62.7, y = -8.9, dx = 0, dy = 0)
p1 <- ggplot() + geom_raster(data = b, aes(x, y, fill = v)) +
  scale_fill_viridis_c(option = "cividis", name = "Depth (m)") +
  geom_sf(data = zones, aes(colour = label, linetype = label), fill = NA, linewidth = 0.7) +
  scale_colour_manual(values = okabe, name = zone_title) + scale_linetype_manual(values = ltys, name = zone_title) +
  guides(colour = guide_legend(title.position = "top", nrow = 2), linetype = guide_legend(title.position = "top", nrow = 2)) +
  leader(i1) + marker(i1) + lab_box(i1, 2.9) + lab_box(b1, 2.7, "italic") + lab_box(jma1, 2.9, "bold", "#D55E00") +
  deg_axes + coord_sf(crs = 4326, expand = FALSE) + labs(x = "Longitude", y = "Latitude") + theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical", panel.background = element_rect(fill = "grey90", colour = NA))
ggsave(file.path(OUT_FIG, "Fig1_study_area.png"), p1, width = 8, height = 8.8, dpi = 300, bg = "white", device = dev)
# Figure 7, predicted suitability from Section 7 for (a) the Chagos Archipelago and (b) the rest of the EEZ and the JMA
pr <- rast(file.path(OUT_R, "v2_pred_full_resolution.tif"))
cut_to <- function(z) { v <- vect(st_union(z)); mask(crop(pr, v), v) }
to_df3 <- function(r) { d <- as.data.frame(aggregate(r, 3, "mean", na.rm = TRUE), xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"; d }
mapp <- function(r, zsub, tag) ggplot() + geom_raster(data = to_df3(r), aes(x, y, fill = v)) +
  geom_sf(data = zsub, aes(linetype = label), fill = NA, colour = "grey30", linewidth = 0.4, show.legend = FALSE) +
  scale_linetype_manual(values = ltys) + scale_fill_viridis_c(name = "Predicted\nsuitability", limits = c(0, 1)) +
  deg_axes + coord_sf(crs = 4326, expand = FALSE) + labs(x = "Longitude", y = "Latitude", tag = tag) + theme_minimal(base_size = 10) +
  theme(panel.background = element_rect(fill = "grey92", colour = NA))
za <- zones[zones$name == "Chagos", ]; zb <- zones[zones$name != "Chagos", ]
ia <- nudge(isl[isl$name %in% c("Diego Garcia", "Peros Banhos", "Salomon Islands", "Egmont Islands"), ],
            dx = c(0.9, -0.9, 0.95, -0.95), dy = c(-0.35, 0.3, 0.3, -0.3))
ba <- nudge(banks[banks$name == "Great Chagos Bank", ], 0, 0)
ib <- nudge(isl[isl$name %in% c("Mauritius", "R\u00e9union", "Rodrigues", "Agal\u00e9ga", "St Brandon", "Tromelin"), ],
            dx = c(1.3, -1.1, 0, 1.2, -1.3, 0), dy = c(-0.6, -0.6, 0.7, 0.3, 0.3, 0.7))
bb <- nudge(banks[banks$name != "Great Chagos Bank", ], dx = c(1.9, 1.6), dy = c(0.3, 0))
jb <- data.frame(name = "Joint Management\nArea", x = 63.0, y = -8.3, dx = 0, dy = 0)
pa <- mapp(cut_to(za), za, "a") + leader(ia) + marker(ia) + lab_box(ia, 2.6) + lab_box(ba, 2.5, "italic")
pb <- mapp(cut_to(zb), zb, "b") + leader(ib) + marker(ib) + lab_box(ib, 2.6) + leader(bb, 0.55) +
  lab_box(bb, 2.5, "italic") + lab_box(jb, 2.6, "bold", "#D55E00")
ggsave(file.path(OUT_FIG, "Fig7_suitability_maps.png"), pa / pb + plot_layout(guides = "collect", heights = c(1, 1.1)),
       width = 8, height = 12, dpi = 300, bg = "white", device = dev)
