#' Regional (gridded) combined index with uncertainty, from the bootstrap replicates
#'
#' The uncertainty counterpart of `computeRegionalIndex()` / `computeGriddedCombinedIndex()`, on the RAW grid (no smoothing: the
#' smoothed map stays a display of the point estimate, made from the final raw result). Inputs are the per-replicate cell means that
#' models_Monitor writes (`<species>/regional/regional_means_<km>km.rds`: mean probability of the German 200 m pixels of each coarse
#' cell, per year and replicate). Inside EVERY replicate, exactly as the baseline does it:
#' - species index per cell = 100 x cell mean in year t / cell mean in the baseline year; a species is left out of a cell where its
#'   baseline mean is below `minBaseline` (the same floor as `computeGriddedCombinedIndex()`);
#' - combined index per cell = geometric mean of the species indices (`exp(mean(log()))`, species without a value left out).
#' Then, over the replicates (replicate 0 = the main models is NOT part of any interval; it is the parity check against the baseline):
#' percentile summaries per year, change between years, trend per decade, and the shares of replicates decreasing / increasing.
#'
#' Output files in `outputDir` (one set per grid size `<km>`):
#' - `regional_index_<km>km_unc_<year>.tif` 5 layers: mean, sd, lwr, upr, width of the combined index;
#' - `regional_index_<km>km_unc_change_<vsBaseline|vs5YearsAgo|vsLastYear>.tif` 7 layers: deltaMean, deltaSd, deltaLwr, deltaUpr,
#'   deltaWidth, shareDecrease, shareIncrease (index of the current year minus index of the reference year, per replicate);
#' - `regional_index_<km>km_unc_trend_per_decade.tif` 7 layers: slopeMean ... slopeWidth, shareDecrease, shareIncrease;
#' - `regional_index_<km>km_coverage.tif` share of the cell's 200 m pixels inside the German outline (border cells < 1);
#' - `regional_index_<km>km_replicates.rds` the per-replicate index array (cells x years x replicates);
#' - `regional_parity_<km>km.txt` replicate 0 against the baseline `regional_index_<km>km_raw.tif`, if `baselineRegionalDir` is given.
#'
#' @param species Character. Species of the combined index (all must have regional_means files; replicates common to all are used).
#' @param uncertaintyDir models_Monitor's uncertainty folder.
#' @param outputDir Where the files are written (created).
#' @param baselineYear,currentYear Integer. Baseline of the index; current year of the change maps.
#' @param cellSizesM Numeric. Grid sizes in metres (the project setting `sharedRegionalCellSizesM`).
#' @param probs Numeric length 2. Interval percentiles.
#' @param minBaseline Numeric. Per-cell baseline floor; MUST be the same value as in the baseline regional index.
#' @param baselineRegionalDir Folder of the baseline `regional_index_<km>km_raw.tif` (parity check), or NULL.
#' @return Invisibly, a list by grid size of the combined index array (cells x years x replicates).
computeRegionalIndexUncertainty <- function(species, uncertaintyDir, outputDir, baselineYear, currentYear, cellSizesM,
                                             probs = c(0.05, 0.95), minBaseline = 1e-6, baselineRegionalDir = NULL) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  res <- list()
  for (cs in cellSizesM) {
    km <- cs / 1000
    fs <- file.path(uncertaintyDir, gsub(" ", "_", species), "regional", sprintf("regional_means_%dkm.rds", km))
    have <- file.exists(fs)
    if (any(!have)) warning("Regional index ", km, " km: no regional_means for ", paste(species[!have], collapse = ", "), call. = FALSE)
    if (!any(have)) next
    parts <- lapply(fs[have], readRDS)
    ids <- Reduce(intersect, lapply(parts, `[[`, "ids"))
    yrs <- Reduce(intersect, lapply(parts, `[[`, "years"))
    if (!(baselineYear %in% yrs)) { message("Regional index ", km, " km: baseline year ", baselineYear, " not available -- skipped."); next }
    g <- parts[[1]]$grid; cells <- parts[[1]]$cells
    stopifnot(all(vapply(parts, function(p) identical(p$cells, cells), logical(1))))
    comb <- regionalCombineSpecies(lapply(parts, function(p) p$mean[, match(yrs, p$years), match(ids, p$ids), drop = FALSE]),
                                   baseIdx = match(baselineYear, yrs), minBaseline = minBaseline)
    dimnames(comb) <- list(NULL, yrs, ids)
    isMain <- ids == 0; reps3 <- comb[, , !isMain, drop = FALSE]
    message("Regional index ", km, " km: ", length(cells), " cells | ", length(yrs), " years | ", dim(reps3)[3], " replicates | ",
            sum(have), " species")
    rasterOf <- function(mat) regionalRaster(g, cells, mat)
    write <- function(r, name) terra::writeRaster(r, file.path(outputDir, sprintf("regional_index_%dkm_%s.tif", km, name)),
                                                  datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
    for (yr in yrs) write(rasterOf(regionalSummary(reps3[, as.character(yr), ], probs)), sprintf("unc_%d", yr))
    comps <- list(vsBaseline = baselineYear, vs5YearsAgo = currentYear - 5L, vsLastYear = currentYear - 1L)
    if (currentYear %in% yrs) for (nm in names(comps)) if (comps[[nm]] %in% yrs)
      write(rasterOf(regionalChange(reps3[, as.character(currentYear), ] - reps3[, as.character(comps[[nm]]), ], probs, "delta")),
            sprintf("unc_change_%s", nm))
    if (length(yrs) >= 5) {
      w <- (yrs - mean(yrs)) / sum((yrs - mean(yrs))^2) * 10                     # slope per decade
      slope <- matrix(0, dim(reps3)[1], dim(reps3)[3])
      for (i in seq_along(yrs)) slope <- slope + w[i] * reps3[, i, ]
      write(rasterOf(regionalChange(slope, probs, "slope")), "unc_trend_per_decade")
    }
    write(rasterOf(cbind(coverage = parts[[1]]$nGerman / parts[[1]]$fact^2)), "coverage")
    saveRDS(list(cellSizeM = cs, cells = cells, years = yrs, ids = ids, index = comb, grid = g), file.path(outputDir, sprintf("regional_index_%dkm_replicates.rds", km)))
    if (!is.null(baselineRegionalDir) && any(isMain)) regionalParity(comb[, , isMain, drop = TRUE], cells, yrs, g, km,
                                                      file.path(baselineRegionalDir, sprintf("regional_index_%dkm_raw.tif", km)),
                                                      file.path(outputDir, sprintf("regional_parity_%dkm.txt", km)))
    res[[as.character(cs)]] <- comb
  }
  invisible(res)
}

#' Combined index per cell, year and replicate from the per-species cell means (the matrix version of
#' `computeGriddedCombinedIndex()`): species index = 100 x mean / baseline mean (left out where the baseline is below `minBaseline`),
#' combined = geometric mean over species.
#' @param meansBySpecies List of arrays cells x years x replicates.
#' @param baseIdx Index of the baseline year along the year dimension.
#' @return Array cells x years x replicates, NA where no species contributes.
regionalCombineSpecies <- function(meansBySpecies, baseIdx, minBaseline = 1e-6) {
  d <- dim(meansBySpecies[[1]]); sumLog <- array(0, d); n <- array(0L, d)
  for (A in meansBySpecies) {
    base <- A[, baseIdx, , drop = TRUE]; base[!(base >= minBaseline)] <- NA_real_
    for (y in seq_len(d[2])) {
      lg <- log(100 * A[, y, ] / base)
      okv <- !is.na(lg)
      s <- sumLog[, y, ]; s[okv] <- s[okv] + lg[okv]; sumLog[, y, ] <- s
      nn <- n[, y, ]; nn[okv] <- nn[okv] + 1L; n[, y, ] <- nn
    }
  }
  out <- exp(sumLog / n); out[n == 0L] <- NA_real_
  out
}

#' Rows = cells, columns = replicates -> mean, sd, lwr, upr, width over the replicates (NA ignored)
regionalSummary <- function(m, probs = c(0.05, 0.95)) {
  m <- matrix(m, nrow = nrow(m))
  q <- if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowQuantiles(m, probs = probs, na.rm = TRUE, useNames = FALSE)
       else t(apply(m, 1, stats::quantile, probs = probs, na.rm = TRUE, names = FALSE))
  sdv <- apply(m, 1, stats::sd, na.rm = TRUE)
  out <- cbind(mean = rowMeans(m, na.rm = TRUE), sd = sdv, lwr = q[, 1], upr = q[, 2], width = q[, 2] - q[, 1])
  out[is.nan(out)] <- NA_real_
  out
}

#' Same plus the shares of replicates with a decrease / an increase (change or trend per replicate in the columns)
regionalChange <- function(d, probs = c(0.05, 0.95), prefix = "delta") {
  d <- matrix(d, nrow = nrow(d)); s <- regionalSummary(d, probs)
  out <- cbind(s, shareDecrease = rowMeans(d < 0, na.rm = TRUE), shareIncrease = rowMeans(d > 0, na.rm = TRUE))
  out[is.nan(out)] <- NA_real_
  colnames(out) <- c(paste0(prefix, c("Mean", "Sd", "Lwr", "Upr", "Width")), "shareDecrease", "shareIncrease")
  out
}

#' Raster of the stored grid with `mat` (cells x layers) written at the stored cell numbers
regionalRaster <- function(g, cells, mat) {
  r <- terra::rast(nrows = g$nrow, ncols = g$ncol, ext = terra::ext(g$ext), crs = g$crs, nlyrs = ncol(mat))
  v <- matrix(NA_real_, terra::ncell(r), ncol(mat)); v[cells, ] <- mat
  terra::values(r) <- v; names(r) <- colnames(mat)
  r
}

#' Replicate 0 (the main models, through the uncertainty machinery) against the baseline raw regional index
regionalParity <- function(main, cells, yrs, g, km, baselineFile, outFile) {
  if (!file.exists(baselineFile)) { message("Regional parity ", km, " km skipped: ", baselineFile, " not found"); return(invisible(NULL)) }
  b <- terra::rast(baselineFile)
  common <- intersect(as.character(yrs), names(b))
  if (!length(common)) { message("Regional parity ", km, " km skipped: no common years"); return(invisible(NULL)) }
  bv <- terra::values(b[[common]])[cells, , drop = FALSE]; mv <- main[, common, drop = FALSE]
  both <- !is.na(bv) & !is.na(mv); dif <- abs(bv - mv)[both]
  lines <- c(sprintf("Regional index %d km | replicate 0 vs baseline %s | years %s", km, basename(baselineFile), paste(common, collapse = ",")),
             sprintf("  cells x years: valid in both = %d, only baseline = %d, only replicate 0 = %d", sum(both), sum(!is.na(bv) & is.na(mv)), sum(is.na(bv) & !is.na(mv))),
             sprintf("  |difference| in index points: median = %.5f, 99th percentile = %.5f, max = %.5f", stats::median(dif), stats::quantile(dif, 0.99), max(dif)))
  bad <- sum(!is.na(bv) & is.na(mv)) + sum(is.na(bv) & !is.na(mv)) > 0.01 * sum(both) || stats::quantile(dif, 0.99) > 0.5
  lines <- c(lines, if (bad) "  RESULT: MISMATCH -- do not trust the regional intervals until this is explained" else "  RESULT: OK")
  writeLines(lines, outFile); message(paste(lines, collapse = "\n"))
  invisible(!bad)
}
