#' Find every completed `metaModel()` run for one species, across run folders
#'
#' `metamodelLabel()` encodes only the three scale resolutions, never
#' species/run name/timestamp -- so two different attempts at the same
#' species/resolutions land in identically-named `metamodel_<label>/`
#' folders under different parent run folders. This globs across every
#' `<candidatesRoot>/*/metamodel_*/` folder to find all of them.
#'
#' @param species Character. Latin species name.
#' @param candidatesRoot Character. Path containing every run's output
#'   folder (e.g. `"outputs"`) -- one level up from each `outputs/<runName>/`.
#' @return A `data.frame` (one row per candidate run, `evalSDM()`'s own
#'   columns plus `species`/`metaDir`/`runFolder`), or `NULL` if no
#'   candidate exists anywhere.
findMetaCandidates <- function(species, candidatesRoot) {
  spClean <- gsub(" ", "_", species)
  perfFiles <- Sys.glob(file.path(candidatesRoot, "*", "metamodel_*",
                                   paste0(spClean, "_perf_meta.rds")))
  if (length(perfFiles) == 0) return(NULL)
  do.call(rbind, lapply(perfFiles, function(f) {
    perf <- readRDS(f)
    perf$species <- species
    perf$metaDir <- dirname(f)
    perf$runFolder <- basename(dirname(dirname(f)))
    perf
  }))
}
