#' Confidence intervals for the species indices and the combined multi-species indices, from the bootstrap replicates
#'
#' Reads what models_Monitor's uncertainty workflow wrote (`<uncertaintyDir>/<Species>/area_mean_replicates.csv`: the
#' country mean of the predicted probability of occurrence, for every replicate and year) and does, INSIDE each replicate:
#'   species index   I_s,t,b = 100 * areaMean_s,t,b / areaMean_s,baseline,b          (as computeSpeciesIndex())
#'   SBI_t,b         = exp(mean over species of log I_s,t,b)                          (as computeCombinedIndex())
#'   Chain_t,b       = 100 at the first year, then * exp(mean over species of log(I_t+1 / I_t))   (as computeChainIndex())
#' and only then takes the mean, median and percentile interval over replicates. A replicate is the SAME bootstrap draw
#' for every species, so the combined indices carry the model uncertainty of all species together. The interval is
#' model uncertainty given the training data -- not the species-resampling interval of computeCombinedIndex().
#'
#' Also stitches the community band pieces (expected richness, mean change in probability) into maps.
#'
#' The two DDA/PECBMS combined indices that need a standard error per species and year are also given, with the SD of the
#' species index ACROSS REPLICATES as that standard error (the existing annual report feeds them the spatial SE of the
#' map instead -- a different quantity): "Analytical" (Gregory et al. 2005 delta method, as computeCombinedIndexAnalytical())
#' and "MSI" (Soldaat et al. 2017 Monte Carlo, as computeCombinedIndexMSI()). Their point estimate is the geometric mean
#' of the replicate-mean species indices; the interval follows `probs` (normal approximation / Monte Carlo percentiles).
#'
#' Not covered here: the regional gridded index maps (computeRegionalIndex() keeps its own method).
#'
#' @param species Character. Species (Latin names) to include.
#' @param uncertaintyDir Character. models_Monitor's uncertainty folder (outputs/<run>/uncertainty).
#' @param reportDir Character. Where the interval tables go (the annual_report folder, next to species_index.csv).
#' @param baselineYear,currentYear Integer. Index baseline; report year of the community change maps.
#' @param probs Numeric length 2. Interval percentiles (default 5th-95th = 90%).
#' @param nBands Integer. Number of bands the community pieces were written in.
#' @return Invisibly, a list with `species` and `combined` (data.frames), or NULL if the baseline year was not mapped.
computeIndexUncertainty <- function(species, uncertaintyDir, reportDir, baselineYear, currentYear,
                                     probs = c(0.05, 0.95), nBands = 16) {
  dir.create(reportDir, recursive = TRUE, showWarnings = FALSE)

  # 1. community maps: stitch the band pieces
  cdir <- file.path(uncertaintyDir, "community")
  pcs <- function(kind, name) file.path(cdir, "pieces", sprintf("%s_%s_band%02d.tif", kind, name, seq_len(nBands)))
  stitch <- function(files, out) {
    files <- files[file.exists(files)]
    if (!length(files)) return(invisible(NULL))
    v <- terra::vrt(files, file.path(tempdir(), paste0("stitch_", basename(out), ".vrt")), overwrite = TRUE)
    names(v) <- names(terra::rast(files[1]))
    terra::writeRaster(v, out, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
    invisible(out)
  }
  stitch(pcs("SRprob", currentYear), file.path(cdir, sprintf("richness_expected_unc_%d.tif", currentYear)))
  for (comp in c("vsBaseline", "vs5YearsAgo", "vsLastYear")) {
    stitch(pcs("meanDeltaP", comp), file.path(cdir, sprintf("community_meanDeltaP_unc_%s.tif", comp)))
    stitch(pcs("bcDissim", comp), file.path(cdir, sprintf("community_turnoverBC_unc_%s.tif", comp)))   # Bray-Curtis turnover, 5 layers
  }

  # 2. area means per species and replicate
  ams <- lapply(species, function(sp) {
    f <- file.path(uncertaintyDir, gsub(" ", "_", sp), "area_mean_replicates.csv")
    if (file.exists(f)) utils::read.csv(f) else { warning("No area means for ", sp, " at ", f, call. = FALSE); NULL }
  })
  ams <- ams[!vapply(ams, is.null, logical(1))]
  if (!length(ams)) stop("computeIndexUncertainty(): no area_mean_replicates.csv found under ", uncertaintyDir)
  am <- do.call(rbind, ams)
  am <- am[am$replicate != 0, ]                       # replicate 0 = the main models (consistency check), not part of an interval
  spp <- unique(am$species); years <- sort(unique(am$year))
  if (!(baselineYear %in% years)) {
    message("computeIndexUncertainty(): baseline year ", baselineYear, " was not mapped -- index intervals skipped.")
    return(invisible(NULL))
  }

  # 3. species indices per replicate -> interval per species and year
  idx <- do.call(rbind, lapply(spp, function(s) {
    d <- am[am$species == s, ]
    base <- d[d$year == baselineYear, c("replicate", "areaMean")]; names(base)[2] <- "baseMean"
    d <- merge(d, base, by = "replicate"); d$index <- 100 * d$areaMean / d$baseMean; d
  }))
  q <- function(x, p) unname(stats::quantile(x, p))
  speciesTable <- do.call(rbind, lapply(split(idx, list(idx$species, idx$year), drop = TRUE), function(d) data.frame(
    species = d$species[1], year = d$year[1], nReplicates = nrow(d),
    areaMeanMean = mean(d$areaMean), areaMeanLwr = q(d$areaMean, probs[1]), areaMeanUpr = q(d$areaMean, probs[2]),
    indexMean = mean(d$index), indexMedian = stats::median(d$index), indexLwr = q(d$index, probs[1]), indexUpr = q(d$index, probs[2]))))
  speciesTable <- speciesTable[order(speciesTable$species, speciesTable$year), ]
  rownames(speciesTable) <- NULL
  utils::write.csv(speciesTable, file.path(reportDir, "species_index_uncertainty.csv"), row.names = FALSE)

  # 4. combined indices per replicate (only replicates present for every species)
  repsAll <- Reduce(intersect, lapply(spp, function(s) unique(am$replicate[am$species == s])))
  message("computeIndexUncertainty(): ", length(spp), " species x ", length(repsAll), " replicates common to all species")
  rows <- lapply(repsAll, function(r) {
    im <- matrix(NA_real_, length(spp), length(years), dimnames = list(spp, years))
    for (s in spp) {
      d <- idx[idx$species == s & idx$replicate == r, ]
      im[s, ] <- d$index[match(years, d$year)]
    }
    sbi <- apply(im, 2, function(v) { v <- v[!is.na(v) & v > 0]; if (length(v)) exp(mean(log(v))) else NA_real_ })
    lam <- log(im[, -1, drop = FALSE] / im[, -ncol(im), drop = FALSE])
    chain <- c(100, 100 * exp(cumsum(colMeans(lam, na.rm = TRUE))))
    data.frame(replicate = r, year = years, SBI = sbi, Chain = chain)
  })
  rep <- do.call(rbind, rows)
  utils::write.csv(rep, file.path(uncertaintyDir, "combined_index_replicates.csv"), row.names = FALSE)
  combined <- do.call(rbind, lapply(c("SBI", "Chain"), function(m) do.call(rbind, lapply(split(rep, rep$year), function(d) data.frame(
    method = m, year = d$year[1], nReplicates = nrow(d), estimate = mean(d[[m]]), median = stats::median(d[[m]]),
    lwr = q(d[[m]], probs[1]), upr = q(d[[m]], probs[2]))))))
  rownames(combined) <- NULL

  # 5. Analytical and MSI: standard error of each species index = SD across replicates
  indexMean <- tapply(idx$index, list(idx$species, idx$year), mean)
  indexSd <- tapply(idx$index, list(idx$species, idx$year), stats::sd)
  zq <- stats::qnorm(probs)
  set.seed(42)
  seRows <- lapply(colnames(indexMean), function(yr) {
    ix <- indexMean[, yr]; se <- indexSd[, yr]
    ok <- !is.na(ix) & !is.na(se) & ix > 0
    ix <- ix[ok]; se <- se[ok]; nT <- length(ix)
    if (nT == 0) return(NULL)
    est <- exp(mean(log(ix)))
    seA <- sqrt((est^2 / nT^2) * sum(se^2 / ix^2))
    sims <- vapply(seq_len(1000), function(i) exp(mean(stats::rnorm(nT, mean = log(ix), sd = se / ix))), numeric(1))
    rbind(data.frame(method = "Analytical", year = as.integer(yr), nReplicates = length(repsAll), estimate = est, median = NA_real_,
                     lwr = est + zq[1] * seA, upr = est + zq[2] * seA),
          data.frame(method = "MSI", year = as.integer(yr), nReplicates = length(repsAll), estimate = est, median = stats::median(sims),
                     lwr = q(sims, probs[1]), upr = q(sims, probs[2])))
  })
  combined <- rbind(combined, do.call(rbind, seRows))
  combined <- combined[order(match(combined$method, c("SBI", "Chain", "Analytical", "MSI")), combined$year), ]
  rownames(combined) <- NULL
  utils::write.csv(combined, file.path(reportDir, "combined_index_uncertainty.csv"), row.names = FALSE)
  message("Index intervals -> ", file.path(reportDir, "species_index_uncertainty.csv"), " and combined_index_uncertainty.csv")
  invisible(list(species = speciesTable, combined = combined))
}
