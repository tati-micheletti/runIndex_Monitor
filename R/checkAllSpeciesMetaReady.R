#' Check whether every species in a roster has a complete metaModel() output
#'
#' Generalizes `models_Monitor`'s `checkAllScalesReady()` pattern (per-species,
#' per-scale file-existence check) to "per-species, across the WHOLE roster" --
#' `metaModel()`'s own self-trigger only knows when one species' 3 scales are
#' done, not when every OTHER species in a cluster run is also done, which is
#' what the multi-species index/report step actually needs before it can run
#' meaningfully. Used by `runIndex_Monitor`'s `checkAllInputs` polling event
#' (see DECISIONS.md's 2026-09-28 "runIndex_Monitor" entry).
#'
#' @param species Character vector. Latin names of every species expected to
#'   have a complete `metaModel()` output.
#' @param metaDir Character. `metaModel()`'s output directory.
#' @param requiredYears Integer vector. Every year a species must have a
#'   valid `<sp>_meta_suitability_<year>.tif` for.
#' @return A list: `allReady` (logical) and `missing` (character vector of
#'   species with at least one missing/invalid required year, empty if none).
checkAllSpeciesMetaReady <- function(species, metaDir, requiredYears) {
  missing <- Filter(function(sp) {
    spClean <- gsub(" ", "_", sp)
    !all(vapply(requiredYears, function(yr) {
      f <- file.path(metaDir, paste0(spClean, "_meta_suitability_", yr, ".tif"))
      isValidPredictionRaster(f, layerName = "meta_prob")
    }, logical(1)))
  }, species)

  list(allReady = length(missing) == 0, missing = missing)
}
