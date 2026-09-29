## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "runIndex_Monitor",
  description = paste("Builds the multi-species annual report and regional gridded index",
                       "maps from models_Monitor's already-computed metaModel() prediction",
                       "rasters -- no new model fitting happens here. Extracted from the",
                       "former root-level runIndex.R script (see DECISIONS.md's 2026-09-28",
                       "\"runIndex_Monitor\" entry) so it gets SpaDES's own parameter",
                       "validation and an init-time check for whether the requested",
                       "indices/comparisons can actually be computed given which years'",
                       "metaModel() output actually exists."),
  keywords = c("bird monitor", "population index", "annual report", "regional index"),
  authors = structure(list(list(given = "Tati", family = "Micheletti", role = c("aut", "cre"),
                                 email = "tati.micheletti@gmail.com", comment = NULL),
                           list(given = "Lisa", family = "Hildebrand", role = "aut",
                                email = "lisa.hildebrand@ufz.de", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(runIndex_Monitor = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "runIndex_Monitor.Rmd"),
  reqdPkgs = list("PredictiveEcology/SpaDES.core@development (>= 3.2.0)",
                   "terra", "mgcv", "geodata"),
  parameters = bindrows(
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?"),

    ## Species / years ---------------------------------------------------------------
    defineParameter("species", "character", NA_character_, NA, NA,
                    "Latin names of the FULL species roster the checkAllInputs polling",
                    "event waits on (e.g. sharedSpecies) -- no default (errors if unset).",
                    "Every one of these must have a complete metaModel() output before",
                    "the index/report computation actually runs."),
    defineParameter("indexSpecies", "character", NULL, NA, NA,
                    "NULL (default): use `species` -- every polled species is also",
                    "covered by the report/index. Otherwise a subset of `species` to",
                    "actually build the report/index for (e.g. a restricted test run)."),
    defineParameter("baselineYear", "numeric", 2005, NA, NA,
                    "Index baseline (2005 -- the first year with genuinely real",
                    "landscape-scale training data in this pipeline). Used by the",
                    "SBI/Analytical/MSI methods; the Chain method ignores it."),
    defineParameter("allYears", "numeric", NA_real_, NA, NA,
                    "Full year range to build the per-species/combined index time",
                    "series over -- no default (errors if unset); typically",
                    "predictionYears (see models_Monitor's parameter of the same name)."),
    defineParameter("currentYear", "numeric", NA_real_, NA, NA,
                    "The year this annual report is \"for\" -- no default (errors if",
                    "unset); typically max(habitatYears)."),
    defineParameter("restrictedYears", "numeric", NULL, NA, NA,
                    "NULL (default): skip the Chain-index robustness check. Otherwise",
                    "years with real (non-hindcast) ground truth (e.g. habitatYears) --",
                    "also computes the Chain index restricted to just these years."),
    defineParameter("changeThresh", "numeric", 0.05, NA, NA,
                    "Passed to computeChangeMaps() -- prevalence-change threshold."),
    defineParameter("nBoot", "numeric", 999, NA, NA,
                    "Passed to computeCombinedIndex()/computeChainIndex() -- bootstrap reps."),
    defineParameter("nSim", "numeric", 1000, NA, NA,
                    "Passed to computeCombinedIndexMSI() -- Monte Carlo simulations."),
    defineParameter("useBootstrapSE", "logical", FALSE, NA, NA,
                    "Passed to computeAnnualReport() -- see its own docstring for the",
                    "spatial-SE-vs-bootstrap-SE tradeoff."),
    defineParameter("cellSizesM", "numeric", c(10000, 20000, 50000), NA, NA,
                    "Regional gridded-index cell sizes (m). Must match",
                    "sharedRegionalCellSizesM (sharedConfig.R, repo root)."),

    ## Cross-run automatic model selection ---------------------------------------------
    defineParameter("metaCandidatesRoot", "character", NULL, NA, NA,
                    "NULL (default): today's behavior -- one shared metaDir from this",
                    "run's own outputPath(sim), no scanning, no audit file. Otherwise a",
                    "path (e.g. \"outputs\") containing every run's output folder --",
                    "triggers auto-discovery of every completed metaModel() run for",
                    "each species across <metaCandidatesRoot>/*/metamodel_*/ (this run's",
                    "own output is naturally one such folder too), picks the best by",
                    "metaSelectionMetric, writes an audit CSV of every candidate",
                    "compared, and stages the winners before computeIndices runs. See",
                    "selectBestMetaPerSpecies()/buildResolvedMetaDir()."),
    defineParameter("metaSelectionMetric", "character", "D2", NA, NA,
                    "Which evalSDM() column (AUC/TSS/Kappa/Sens/Spec/PCC/D2/thresh)",
                    "decides the winner when metaCandidatesRoot is set. All columns are",
                    "logged in the audit regardless. Default D2 (explained deviance) --",
                    "matches the source paper's own scale-importance methodology and",
                    "metaModel()'s ridge-GLM family."),

    ## Scale resolutions (to reconstruct metamodelLabel(), must match models_Monitor's) --
    defineParameter("climateResolutionM", "numeric", 50000, NA, NA,
                    "Must match models_Monitor's/dataPrep_Monitor's climateResolutionM."),
    defineParameter("habitatResolutionM", "numeric", 200, NA, NA,
                    "Must match models_Monitor's/dataPrep_Monitor's habitatResolutionM."),
    defineParameter("landscapeResolutionM", "numeric", 1000, NA, NA,
                    "Must match models_Monitor's/dataPrep_Monitor's landscapeResolutionM."),

    ## Cluster-run polling -- irrelevant for a normal single-session runMe.R run, where the
    ## very first checkAllInputs check passes immediately (see DECISIONS.md's 2026-09-28
    ## "runIndex_Monitor" entry for why: SpaDES processes same-time(sim) events in
    ## scheduling order, so loadOrder already puts this module's init after every one of
    ## models_Monitor's own events). Only exercised when models_Monitor was run as separate
    ## per-species-per-scale cluster tasks (tools/runClusterTask.R), which finish at
    ## unpredictable wall-clock times relative to each other. --------------------------
    defineParameter("pollIntervalSeconds", "numeric", 300, NA, NA,
                    "Real wall-clock seconds between checkAllInputs rechecks (Sys.sleep()",
                    "-- this pipeline's simulation timeunit is \"year\", already a",
                    "formality per runMe.R's own comment, so a real elapsed-time poll",
                    "can't be expressed as scheduling a future event some years later)."),
    defineParameter("pollTimeoutHours", "numeric", 6, NA, NA,
                    "After this many real wall-clock hours of polling, WARN (never stop)",
                    "which species are still missing a metaModel() output and proceed",
                    "anyway with whatever is available -- avoids polling forever on a",
                    "paid HPC allocation if a cluster task genuinely failed.")
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
  ),
  outputObjects = bindrows(
    createsOutput("annualReport", "list",
                  "computeAnnualReport()'s full return value -- speciesIndex,",
                  "combinedIndex{SBI,Analytical,MSI,Chain[Restricted]}, srMap, changeMaps."),
    createsOutput("regionalIndex", "list",
                  "computeRegionalIndex()'s full return value -- gridded + smoothed",
                  "regional index rasters per cell size in cellSizesM.")
  )
))

doEvent.runIndex_Monitor = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      if (identical(P(sim)$species, NA_character_)) {
        stop("runIndex_Monitor's species parameter must be supplied explicitly ",
             "(e.g. sharedSpecies from sharedConfig.R) -- no default.")
      }
      if (identical(P(sim)$allYears, NA_real_) || identical(P(sim)$currentYear, NA_real_)) {
        stop("runIndex_Monitor's allYears/currentYear parameters must be supplied ",
             "explicitly (e.g. predictionYears/max(habitatYears)) -- no default.")
      }
      mod$pollStartTime <- Sys.time()
      sim <- scheduleEvent(sim, time(sim), "runIndex_Monitor", "checkAllInputs")
    },

    checkAllInputs = {
      # ! ----- EDIT BELOW ----- ! #
      resolutionsM <- c(europe = P(sim)$climateResolutionM,
                         habitat = P(sim)$habitatResolutionM,
                         landscape = P(sim)$landscapeResolutionM)
      metaDir <- file.path(outputPath(sim), metamodelLabel(resolutionsM))

      # Whichever years the report/index actually needs real metaModel()
      # output for: currentYear + restrictedYears (allYears' broader series
      # tolerates gaps as NA -- see validateIndexYears() -- so isn't a
      # blocking readiness condition the way these are).
      requiredYears <- sort(unique(c(P(sim)$currentYear, P(sim)$restrictedYears)))
      ready <- checkAllSpeciesMetaReady(P(sim)$species, metaDir, requiredYears)

      if (!ready$allReady) {
        elapsedHours <- as.numeric(difftime(Sys.time(), mod$pollStartTime, units = "hours"))
        if (elapsedHours >= P(sim)$pollTimeoutHours) {
          warning("runIndex_Monitor: timed out after ", P(sim)$pollTimeoutHours,
                  "h waiting for metaModel() output -- proceeding with whatever ",
                  "species are ready. Still missing: ",
                  paste(ready$missing, collapse = ", "), call. = FALSE)
        } else {
          message("runIndex_Monitor: ", length(P(sim)$species) - length(ready$missing), "/",
                  length(P(sim)$species), " species ready -- rechecking in ",
                  P(sim)$pollIntervalSeconds, "s. Still missing: ",
                  paste(ready$missing, collapse = ", "))
          Sys.sleep(P(sim)$pollIntervalSeconds)
          sim <- scheduleEvent(sim, time(sim), "runIndex_Monitor", "checkAllInputs")
          return(invisible(sim))
        }
      }

      sim <- scheduleEvent(sim, time(sim), "runIndex_Monitor", "computeIndices")
      # ! ----- STOP EDITING ----- ! #
    },

    computeIndices = {
      # ! ----- EDIT BELOW ----- ! #
      resolutionsM <- c(europe = P(sim)$climateResolutionM,
                         habitat = P(sim)$habitatResolutionM,
                         landscape = P(sim)$landscapeResolutionM)
      metaDir <- file.path(outputPath(sim), metamodelLabel(resolutionsM))
      if (!dir.exists(metaDir)) {
        stop("runIndex_Monitor: no metaModel() output found at ", metaDir,
             " -- has models_Monitor's metaModel event run in this session?")
      }

      indexSpecies <- if (is.null(P(sim)$indexSpecies)) P(sim)$species else P(sim)$indexSpecies

      effectiveMetaDir <- metaDir
      if (!is.null(P(sim)$metaCandidatesRoot)) {
        selection <- selectBestMetaPerSpecies(indexSpecies, P(sim)$metaCandidatesRoot,
                                               P(sim)$metaSelectionMetric)
        effectiveMetaDir <- file.path(outputPath(sim), "metamodel_resolved")
        dir.create(effectiveMetaDir, recursive = TRUE, showWarnings = FALSE)
        auditPath <- file.path(effectiveMetaDir, "meta_model_selection_audit.csv")
        write.csv(selection$audit, auditPath, row.names = FALSE)
        buildResolvedMetaDir(selection$winners, effectiveMetaDir)
        message("runIndex_Monitor: auto-selected best metaModel() per species by ",
                P(sim)$metaSelectionMetric, " -- audit: ", auditPath)
      }

      validateIndexYears(species = indexSpecies, allYears = P(sim)$allYears,
                          baselineYear = P(sim)$baselineYear, currentYear = P(sim)$currentYear,
                          restrictedYears = P(sim)$restrictedYears, metaDir = effectiveMetaDir)

      reportDir <- file.path(outputPath(sim), "annual_report")
      regionalDir <- file.path(outputPath(sim), "regional_index")

      # Same GADM level-0 fetch/cache pattern as bootstrapMetaModelTrend()'s
      # (models_Monitor) -- avoids cells straddling the border silently
      # averaging in non-German source data.
      germanyBoundary <- NULL
      if (requireNamespace("geodata", quietly = TRUE)) {
        gadmCacheDir <- file.path(inputPath(sim), "predictors", "raw", "gadm")
        dir.create(gadmCacheDir, recursive = TRUE, showWarnings = FALSE)
        germanyBoundary <- geodata::gadm(country = "DEU", level = 0, path = gadmCacheDir)
      } else {
        warning("'geodata' package not installed -- computeRegionalIndex() will run ",
                "without a country boundary (cells near the border may average in ",
                "non-German source data).", call. = FALSE)
      }

      message("=== Annual report (", P(sim)$currentYear, ") ===")
      sim$annualReport <- computeAnnualReport(
        species = indexSpecies,
        baselineYear = P(sim)$baselineYear,
        currentYear = P(sim)$currentYear,
        allYears = P(sim)$allYears,
        metaDir = effectiveMetaDir,
        outputDir = reportDir,
        restrictedYears = P(sim)$restrictedYears,
        changeThresh = P(sim)$changeThresh,
        nBoot = P(sim)$nBoot,
        nSim = P(sim)$nSim,
        useBootstrapSE = P(sim)$useBootstrapSE)

      message("\n=== Regional index maps ===")
      sim$regionalIndex <- computeRegionalIndex(
        species = indexSpecies,
        years = P(sim)$allYears,
        baselineYear = P(sim)$baselineYear,
        metaDir = effectiveMetaDir,
        outputDir = regionalDir,
        cellSizesM = P(sim)$cellSizesM,
        countryBoundary = germanyBoundary)

      message("\nDone. Report -> ", reportDir, " | Regional maps -> ", regionalDir)
      # ! ----- STOP EDITING ----- ! #
    },

    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
