#' Stage each species' winning `metaModel()` output into one directory
#'
#' Copies each species' `<spClean>_meta_suitability_*.tif` files (every
#' year present) and `<spClean>_meta_trend_boot.rds` (if present) from its
#' winning run's `metaDir` into `stagingDir`, so `validateIndexYears()`/
#' `computeAnnualReport()`/`computeRegionalIndex()` can read one shared
#' directory exactly as they do today, regardless of how many different
#' runs the winners actually came from.
#'
#' @param winners Named list, species -> that species' winning `metaDir`
#'   (as returned by `selectBestMetaPerSpecies()`'s own `winners` element).
#' @param stagingDir Character. Directory to copy into (created if needed).
#' @return Invisibly, `stagingDir`.
buildResolvedMetaDir <- function(winners, stagingDir) {
  dir.create(stagingDir, recursive = TRUE, showWarnings = FALSE)
  for (sp in names(winners)) {
    spClean <- gsub(" ", "_", sp)
    pattern <- paste0("^", spClean, "_meta_(suitability_.*\\.tif|trend_boot\\.rds)$")
    files <- list.files(winners[[sp]], pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
      warning("buildResolvedMetaDir: no meta output files found for ", sp,
               " under ", winners[[sp]], call. = FALSE)
      next
    }
    file.copy(files, stagingDir, overwrite = TRUE)
  }
  invisible(stagingDir)
}
