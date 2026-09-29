#' Pick the best-performing `metaModel()` run per species, with a full audit trail
#'
#' Compares every candidate run `findMetaCandidates()` finds for each
#' species and picks the one maximizing `selectionMetric`
#' (`evalSDM()`'s own column names: `AUC`/`TSS`/`Kappa`/`Sens`/`Spec`/
#' `PCC`/`D2`/`thresh`). Purely performance-metric-based -- does NOT also
#' weigh year-completeness of each candidate (that's `validateIndexYears()`'s
#' job, run afterward against whichever run gets selected here).
#'
#' @param species Character vector. Species to select a run for.
#' @param candidatesRoot Character. Passed to `findMetaCandidates()`.
#' @param selectionMetric Character, default `"D2"` (explained deviance --
#'   matches the source paper's own scale-importance methodology, and
#'   `metaModel()` is a ridge GLM so D2 reflects overall calibration/fit,
#'   not just discrimination).
#' @return A list: `audit` (data.frame, every candidate for every species,
#'   with a `selected` logical column) and `winners` (named list,
#'   species -> that species' winning `metaDir`).
selectBestMetaPerSpecies <- function(species, candidatesRoot, selectionMetric = "D2") {
  perSpecies <- lapply(species, findMetaCandidates, candidatesRoot = candidatesRoot)
  names(perSpecies) <- species
  found <- perSpecies[!vapply(perSpecies, is.null, logical(1))]
  missing <- setdiff(species, names(found))
  if (length(missing) > 0) {
    warning("No metaModel() candidates found anywhere under ", candidatesRoot,
            " for: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(found) == 0) {
    return(list(audit = data.frame(), winners = list()))
  }

  audit <- do.call(rbind, found)
  audit$selected <- FALSE
  for (sp in names(found)) {
    idx <- which(audit$species == sp)
    best <- idx[which.max(audit[[selectionMetric]][idx])]
    audit$selected[best] <- TRUE
  }
  winners <- stats::setNames(audit$metaDir[audit$selected], audit$species[audit$selected])
  list(audit = audit, winners = as.list(winners))
}
