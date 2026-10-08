#' Prevalence-change and gain/loss maps between two years
#'
#' Per Zbinden et al. (2005)'s own caution -- a combined/community index
#' can mask species-specific heterogeneity (one species' strong increase
#' can offset another's strong decline), so per-species indices must
#' always be examined alongside it, never the combined index in
#' isolation -- this returns both per-species maps and their community
#' aggregate, for every comparison.
#'
#' Three map types, all cell-by-cell raster algebra on existing
#' `metaModel()` output -- no new model fitting:
#' - **Prevalence change** (continuous, no threshold): `deltaP = P_current - P_ref`.
#' - **Gain/loss** (categorical, from the already-thresholded `binary`
#'   band): 0 = stable absent, 1 = loss, 2 = gain, 3 = stable present.
#' - **Stable/increase/decrease** (categorical, from `deltaP` and
#'   `changeThresh`): -1 = decrease, 0 = stable, 1 = increase.
#'
#' @param species Character vector of Latin species names.
#' @param yearRef Integer. The earlier ("reference") year.
#' @param yearCurrent Integer. The later ("current") year.
#' @param metaDir Character. Directory holding metaModel()'s output.
#' @param changeThresh Numeric. Minimum absolute `deltaP` to call a cell
#'   "increase"/"decrease" rather than "stable". Default 0.05 -- a
#'   judgment call, state it explicitly in any methods write-up.
#' @return List with:
#'   - `perSpecies`: named list (by species) of 3-layer SpatRasters
#'     (`deltaP`, `gainLoss`, `stableIncDec`), or `NULL` for a species
#'     missing either year.
#'   - `community`: a 5-layer SpatRaster (`meanDeltaP`, `netGainLoss`,
#'     `gainCount`, `lossCount`, `turnoverBC`), aggregated across all species
#'     with valid data for both years. `turnoverBC` is the Bray-Curtis
#'     dissimilarity of the expected communities of the two years,
#'     sum(|deltaP|) / sum(P_ref + P_current): 0 = unchanged, 1 = no species
#'     in common. Unlike `meanDeltaP` (gains and losses cancel) and expected
#'     richness (a swap of species leaves the sum unchanged), it shows turnover.
computeChangeMaps <- function(species, yearRef, yearCurrent, metaDir, changeThresh = 0.05) {
  perSpecies <- list()
  deltaPLayers <- list()
  absLayers <- list()     # |deltaP| (Bray-Curtis numerator)
  totalLayers <- list()   # P_ref + P_current (Bray-Curtis denominator)
  signedLayers <- list()
  gainLayers <- list()
  lossLayers <- list()

  readBands <- function(sp, yr) {
    spClean <- gsub(" ", "_", sp)
    f <- file.path(metaDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))
    if (!file.exists(f)) return(NULL)
    tryCatch(terra::rast(f), error = function(e) NULL)
  }

  for (sp in species) {
    rRef <- readBands(sp, yearRef)
    rCur <- readBands(sp, yearCurrent)

    if (is.null(rRef) || is.null(rCur) ||
        !all(c("meta_prob", "binary") %in% names(rRef)) ||
        !all(c("meta_prob", "binary") %in% names(rCur))) {
      perSpecies[[sp]] <- NULL
      next
    }

    # Every computed (not read-from-disk) raster below is force-materialized
    # fully into memory immediately (terra::values(x) <- terra::values(x))
    # -- a raster arithmetic result otherwise stays lazily linked to its
    # source SpatRasters (rRef/rCur here), which can silently go invalid
    # later in a long-lived session with many other raster reads/writes in
    # between (found via a real failure chaining this function into the
    # same SpaDES session as models_Monitor; it worked fine called
    # standalone in the fresh, short-lived R session the pre-existing
    # root-level runIndex.R script always ran it in -- see DECISIONS.md's
    # 2026-09-28 "runIndex_Monitor" entry).
    deltaP <- rCur[["meta_prob"]] - rRef[["meta_prob"]]
    names(deltaP) <- "deltaP"
    terra::values(deltaP) <- terra::values(deltaP)

    gainLoss <- rRef[["binary"]] + 2 * rCur[["binary"]]
    names(gainLoss) <- "gainLoss"
    terra::values(gainLoss) <- terra::values(gainLoss)

    stableIncDec <- terra::classify(
      deltaP,
      rcl = matrix(c(-Inf, -changeThresh, -1,
                      -changeThresh, changeThresh, 0,
                      changeThresh, Inf, 1),
                    ncol = 3, byrow = TRUE)
    )
    names(stableIncDec) <- "stableIncDec"
    terra::values(stableIncDec) <- terra::values(stableIncDec)

    perSpecies[[sp]] <- combineLayersSafely(list(deltaP, gainLoss, stableIncDec))

    deltaPLayers[[sp]] <- deltaP
    absP <- abs(deltaP); terra::values(absP) <- terra::values(absP)
    absLayers[[sp]] <- absP
    totP <- rRef[["meta_prob"]] + rCur[["meta_prob"]]; terra::values(totP) <- terra::values(totP)
    totalLayers[[sp]] <- totP
    signedGL <- (gainLoss == 2) - (gainLoss == 1)  # +1 gain, -1 loss, 0 else
    terra::values(signedGL) <- terra::values(signedGL)
    signedLayers[[sp]] <- signedGL
    gainGL <- (gainLoss == 2); terra::values(gainGL) <- terra::values(gainGL)
    gainLayers[[sp]] <- gainGL
    lossGL <- (gainLoss == 1); terra::values(lossGL) <- terra::values(lossGL)
    lossLayers[[sp]] <- lossGL
  }

  if (length(deltaPLayers) == 0) {
    warning("No species had valid data for both ", yearRef, " and ", yearCurrent,
            " -- community layer is NULL")
    community <- NULL
  } else {
    # terra::mean()/terra::app() called explicitly, namespaced -- NOT the
    # bare mean()/sum() generics. Found via a real failure chaining this
    # function into the same SpaDES session as models_Monitor: something
    # about that session's package loading leaves the base mean()/sum()
    # generics' S4 dispatch to terra's SpatRaster methods broken (bare
    # mean(spatRaster) silently falls through to mean.default() and
    # returns a numeric NA instead of a raster) -- terra's OWN explicitly-
    # namespaced functions sidestep that dispatch entirely, regardless of
    # its root cause. See DECISIONS.md's 2026-09-28 "runIndex_Monitor" entry.
    meanDeltaP <- terra::mean(combineLayersSafely(deltaPLayers), na.rm = TRUE)
    names(meanDeltaP) <- "meanDeltaP"
    terra::values(meanDeltaP) <- terra::values(meanDeltaP)
    netGainLoss <- terra::app(combineLayersSafely(signedLayers), fun = "sum", na.rm = TRUE)
    names(netGainLoss) <- "netGainLoss"
    terra::values(netGainLoss) <- terra::values(netGainLoss)
    gainCount <- terra::app(combineLayersSafely(gainLayers), fun = "sum", na.rm = TRUE)
    names(gainCount) <- "gainCount"
    terra::values(gainCount) <- terra::values(gainCount)
    lossCount <- terra::app(combineLayersSafely(lossLayers), fun = "sum", na.rm = TRUE)
    names(lossCount) <- "lossCount"
    terra::values(lossCount) <- terra::values(lossCount)
    sumAbs <- terra::app(combineLayersSafely(absLayers), fun = "sum", na.rm = TRUE)
    sumTot <- terra::app(combineLayersSafely(totalLayers), fun = "sum", na.rm = TRUE)
    turnoverBC <- sumAbs / terra::ifel(sumTot < 1e-12, 1e-12, sumTot)
    names(turnoverBC) <- "turnoverBC"
    terra::values(turnoverBC) <- terra::values(turnoverBC)
    community <- combineLayersSafely(list(meanDeltaP, netGainLoss, gainCount, lossCount, turnoverBC))
  }

  list(perSpecies = perSpecies, community = community)
}
