#' Warn (never stop) about years a requested index/report can't cover
#'
#' `computeChangeMaps()`'s `vsBaseline`/`vs5YearsAgo`/`vsLastYear`
#' comparisons already degrade gracefully to `NULL`/`NA` when a needed
#' year's `metaModel()` output is missing -- silently, discovered only by
#' inspecting the output. This surfaces the same gaps up front, by species,
#' naming exactly which year and which comparison is affected, before the
#' (possibly long-running) computation starts. Always `warning()`s, never
#' `stop()`s -- a sparse or restricted `predictionYears`/`allYears` is a
#' deliberate, valid choice (e.g. a caching/scale test), not an error.
#'
#' @param species Character vector. Latin names to check.
#' @param allYears Integer vector. Full year range the per-species/combined
#'   index series is built over.
#' @param baselineYear Integer. Index baseline year.
#' @param currentYear Integer. The year this report is "for".
#' @param restrictedYears Integer vector or NULL. Years for the Chain-index
#'   robustness check.
#' @param metaDir Character. `metaModel()`'s output directory.
#' @return Invisible NULL. Called for its `warning()` side effects.
validateIndexYears <- function(species, allYears, baselineYear, currentYear,
                                restrictedYears, metaDir) {
  namedYears <- list(baselineYear = baselineYear, currentYear = currentYear,
                      vsLastYear = currentYear - 1, vs5YearsAgo = currentYear - 5)

  for (sp in species) {
    spClean <- gsub(" ", "_", sp)
    hasYear <- function(yr) {
      isValidPredictionRaster(file.path(metaDir, paste0(spClean, "_meta_suitability_", yr, ".tif")),
                               layerName = "meta_prob")
    }

    missingNamed <- Filter(function(nm) !hasYear(namedYears[[nm]]), names(namedYears))
    missingSeries <- Filter(function(yr) !hasYear(yr), allYears)
    missingRestricted <- if (!is.null(restrictedYears)) {
      Filter(function(yr) !hasYear(yr), restrictedYears)
    } else {
      integer(0)
    }

    msgParts <- character(0)
    for (nm in missingNamed) {
      msgParts <- c(msgParts, sprintf(
        "%s needs year %d's metaModel() output, which is missing -- that comparison will return NA.",
        nm, namedYears[[nm]]))
    }
    if (length(missingSeries) > 0) {
      msgParts <- c(msgParts, sprintf(
        "the per-species/combined index series is missing %d of %d requested years (%s) -- those points will be NA.",
        length(missingSeries), length(allYears), paste(missingSeries, collapse = ", ")))
    }
    if (length(missingRestricted) > 0) {
      msgParts <- c(msgParts, sprintf(
        "the Chain-index robustness check is missing %d of %d restrictedYears (%s).",
        length(missingRestricted), length(restrictedYears), paste(missingRestricted, collapse = ", ")))
    }

    if (length(msgParts) > 0) {
      warning(sp, ": ", paste(msgParts, collapse = " "), call. = FALSE)
    }
  }

  invisible(NULL)
}
