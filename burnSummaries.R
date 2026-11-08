defineModule(sim, list(
  name = "burnSummaries",
  description = "Extract simulated fire sizes to generate fire regime statistics",
  keywords = "",
  authors = c(
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre"),
           comment = c(ORCID = "0000-0001-7146-8135"))
  ),
  childModules = character(0),
  version = list(burnSummaries = "1.1.0.9004"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "burnSummaries.Rmd"),
  loadOrder = list(after = c("fireSense_burn", "LandMine", "scfmSpread")),
  ## fireregimetools >= 0.1.0.9003: fetch_nbac_polys()/fetch_nfdb_polys() acquire the national fire
  ## archives with a raised download timeout and verify each extraction against the archive
  ## manifest's sizes instead of merely checking existence (a truncated extraction used to be
  ## silently reused forever). Earlier versions lack the fetchers. >= 0.1.0.9008: they read only
  ## the study area's extent of each national file (FOR-CAST/fireregimetools#2), with that extent
  ## built from the study area's outline so records near a curved or reprojected edge are kept. `archive` is not called here, but
  ## fireregimetools uses libarchive for extraction when it is installed, which handles the zip64
  ## archives R's internal unzip cannot -- so keep it available.
  reqdPkgs = list("archive", "data.table", "dplyr", "httr2", "PredictiveEcology/LandR@development",
                  "FOR-CAST/fireregimetools@main (>= 0.1.0.9008)",
                  "ggplot2", "ggspatial", "googledrive", "patchwork",
                  "reproducible", "sf",
                  ## dirnamesFromSet() (bare-number replicate folders), resolveSimYears(), padYears()
                  "PredictiveEcology/SpaDES.core@development (>= 3.2.1.9057)",
                  "stringr", "terra", "tidyterra"),
  parameters = bindrows(
    defineParameter("dataYear", "numeric", 2020, NA, NA,
                    "data year for the SCANFI stand-age inputs used to seed rstTimeSinceFire (single mode)"),
    defineParameter("fireTimestep", "integer", 1L, NA, NA,
                    "simulation time interval between burn events"),
    defineParameter("historicYears", "integer", c(NA_integer_, NA_integer_), NA, NA,
                    paste("(multi mode) first and last year of the observed fire record used for the historic",
                          "mean annual burn map and the observed row of the fire regime table, e.g.",
                          "`c(1985L, 2024L)`. `NA` (the default) uses every year the record holds.")),
    defineParameter("mode", "character", "single", NA, NA,
                    paste("use 'single' to run part of a simulation;",
                          "use 'multi' to run as part of postprocessing multiple runs.")),
    defineParameter("nbacStartYear", "integer", 1972L, NA, NA,
                    paste("first year of the historic fire record (multi mode) taken from NBAC",
                          "perimeters; earlier years come from NFDB polygons. NBAC begins in 1972.")),
    defineParameter("reps", "integer", 1L:10L, 1L, NA_integer_,
                    paste("in single mode, should be length 1 indicating the replicate id;",
                          "in multi mode, a vector of replicate IDs.")),
    defineParameter("simOutputPath", "character", outputPath(sim), NA, NA,
                    "Directory specifying the location of the simulation outputs."),
    defineParameter("simTimes", "numeric", c(NA, NA), NA, NA,
                    "Simulation start and end times when running in 'multi' mode."),
    defineParameter("summaryInterval", "integer", 100L, NA, NA,
                    "simulation time interval at which to take 'snapshots' used for summary analyses"),
    defineParameter("summaryPeriod", "integer", c(NA_integer_, NA_integer_), NA, NA,
                    paste("lower and upper end of the range of simulation times used for summary analyses.",
                          "`NA` (the default) means the last 30% of the run, which is years 700 to 1000 of a",
                          "1000-year run; a run shorter than 1000 years gets a message with the period chosen.")),
    defineParameter(".plots", "character", "png", NA, NA,
                    "Used by `Plots` function, which can be optionally used here."),
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
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("burnMap", "SpatRaster",
                 desc = paste("Cumulative burn map.", "Required in single mode.")),
    expectsInput("burnSummary", "data.table",
                 desc = paste("Fire summary table from `fireSense_burn` or `scfm`.",
                              "One of `burnSummary` or `fireSizes` is required in single mode.")),
    expectsInput("firePolys", "list", sourceURL = NA,
                 desc = paste("Optional, multi mode. List of `SpatVector` polygons of observed fires,",
                              "each with a `YEAR` column. If not supplied, NBAC perimeters",
                              "(with NFDB backfill) are downloaded.")),
    expectsInput("fireReturnInterval", "SpatRaster",
                 desc = paste("Fire return interval (years), on the grid of the saved `flammableMap`;",
                              "LandMine's input of that name. Multi mode, optional: both burn maps are",
                              "restricted to cells with an interval (not NA, above 0), the cells where",
                              "LandMine ignites fires.")),
    expectsInput("fireSizes", "list",
                 desc = paste("Fire sizes summary tables from LandMine.",
                              "One of `burnSummary` or `fireSizes` is required in single mode.")),
    expectsInput("flammableMap", "SpatRaster",
                 desc = paste("Binary flammability map.",
                              "Required in single mode.")),
    expectsInput("outputsDF",  "data.table",
                 desc = paste("Optional, multi mode. The rbindlisted `outputs(sim)` of all the sims",
                              "being summarised, i.e., all the files that may exist.",
                              "If not supplied, files are read from `outputPath(sim)/rep<NN>`.")),
    expectsInput("nonForest_timeSinceDisturbance",  "SpatRaster",
                 desc = paste("map of time since last burn, with non-flammable pixels receiving `NA`.",
                              "One of `rstTimeSinceFire` or `nonForest_timeSinceDisturbance` is required in single mode.")),
    expectsInput("rstCurrentBurn", "SpatRaster",
                 desc = "Binary raster of fires, 1 meaning 'burned', 0 or NA is non-burned"),
    expectsInput("rstTimeSinceFire", "SpatRaster",
                 desc = paste("map of time since last burn, with non-flammable pixels receiving `NA`.",
                              "One of `rstTimeSinceFire` or `nonForest_timeSinceDisturbance` is required in single mode.")),
    expectsInput("studyArea", "SpatVector",
                 desc = paste("Polygon of the area the fire model burns; LandMine's input of that name.",
                              "Multi mode, optional: both burn maps are restricted to cells whose centre",
                              "is inside it. An `sf` object is also accepted.")),
  ),
  outputObjects = bindrows(
    createsOutput("fireSizes", "data.table",
                  desc = "summary fire sizes table"),
    createsOutput("rstTimeSinceFire", "SpatRaster",
                  desc = "map of time since last burn, with non-flammable pixels receiving `NA`.")
  )
))

## event types
#   - type `init` is required for initialization

doEvent.burnSummaries = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      params(sim)$burnSummaries$summaryPeriod <- resolveSummaryPeriod(P(sim)$summaryPeriod, start(sim), end(sim))

      ## summary output times = seq over the summary period (inlined so this generic
      ## module does not depend on LandWebUtils::analysesOutputsTimes -- which is only
      ## loaded when a LandWeb module like NRV_summary is co-run, not in standalone
      ## mode="multi").
      mod$analysesOutputsTimes <- summaryTimes(P(sim)$summaryPeriod, P(sim)$summaryInterval)

      if (P(sim)$mode == "single") {
        sim <- InitSingle(sim)

        sim <- scheduleEvent(sim, start(sim), "burnSummaries", "update_tsf")
        sim <- scheduleEvent(sim, end(sim), "burnSummaries", "create_fireSizes", .last())

        ## Year 0 is the map the run starts from, so it is saved before update_tsf (priority 5) adds
        ## fireTimestep to it, and after every module's init (SpaDES.core schedules those at
        ## .first()). At .last(), the saved year 0 was the starting map + fireTimestep.
        sim <- scheduleEvent(sim, start(sim), "burnSummaries", "save_single", .first() + 1)

        ## The other saves come at the end of their year. A time equal to start(sim) is left out: a
        ## second year-0 save would overwrite the one above; the summary period start and the end
        ## are each scheduled once, even when they coincide. When the summary period starts at
        ## start(sim), the save above schedules the next snapshot.
        later <- c(P(sim)$summaryPeriod[1], end(sim))
        for (saveTime in unique(later[later != start(sim)])) {
          sim <- scheduleEvent(sim, saveTime, "burnSummaries", "save_single", .last())
        }
      } else if (P(sim)$mode == "multi") {
        sim <- InitMulti(sim)

        ## schedule future event(s)
        sim <- scheduleEvent(sim, start(sim), "burnSummaries", "summary")
        sim <- scheduleEvent(sim, start(sim), "burnSummaries", "plot")
      }
    },
    update_tsf = {
      sim$rstTimeSinceFire[] <- as.integer(sim$rstTimeSinceFire[]) + as.integer(P(sim)$fireTimestep) ## preserves NAs
      sim$rstTimeSinceFire[which(sim$rstCurrentBurn[] == 1)] <- 0L

      ## schedule next event
      sim <- scheduleEvent(sim, time(sim) + P(sim)$fireTimestep, "burnSummaries", "update_tsf")
    },
    create_fireSizes = {
      ## try to guess the replicate being run
      if (length(P(sim)$reps) == 1) {
        repID <- P(sim)$reps
      } else if (grepl("rep", outputPath(sim))) {
        repID <- stringr::str_extract(outputPath(sim), "rep[0-9].*") |>
          stringr::str_remove("rep") |>
          as.integer()
      } else {
        repID <- NA_integer_
      }
      if (!is.null(sim[["burnSummary"]])) {
        fs <- copy(sim[["burnSummary"]])
        fs[, `:=`(
          simArea = P(sim)$.studyAreaName,
          size = N, ## use number of pixels (N) instead of areaBurned
          maxSize = NA_integer_, ## NOTE: no "target" or "expected" fire size w/ scfm nor fS
          rep = repID,
          ## LandMine emits these; scfm and fireSense have no notion of a retried fire. Declared
          ## NA rather than omitted so every fire model writes the SAME schema -- arrow unifies
          ## schemas across the per-replicate parquet files, and a mixed set would fail to open.
          fireID = NA_integer_,
          attempt = NA_integer_,
          targetSize = NA_real_
        )]
        set(fs, NULL, c("areaBurned", "igLoc", "grp", "PolyID"), NULL)
        setcolorder(fs, c("simArea", "rep", "year", "size", "maxSize"))
        setnames(fs, old = c("size", "maxSize"), new = c("simSize", "expSize"))
      } else if (!is.null(sim[["fireSizes"]])) {
        fs <- rbindlist(sim[["fireSizes"]], idcol = "year")
        fs[, `:=`(simArea = P(sim)$.studyAreaName, rep = repID)]
        setcolorder(fs, c("simArea", "rep", "year", "size", "maxSize"))
        setnames(fs, old = c("size", "maxSize"), new = c("simSize", "expSize"))
      } else {
        stop("cumulative fire summary object not found")
      }

      ffs <- file.path(outputPath(sim), "burnSummaries_fireSizes.csv")
      data.table::fwrite(fs, file = ffs)

      sim <- registerOutputs(ffs, sim)

      ## also publish as a parquet partition so multi mode reads all reps as one lazy Arrow dataset
      ## (fireregimetools::open_burn_dataset) instead of rbind-ing per-rep CSVs into memory.
      fireregimetools::write_burn_parquet(
        as.data.frame(fs),
        file.path(outputPath(sim), "burnSummaries_fireSizes"),
        replicate = repID
      )
    },
    summary = {
      sim <- FireSummaries(sim)
    },
    plot = {
      plotFun(sim)
    },
    save_single = {
      padYear <- paddedFloatToChar(time(sim), padL = ceiling(log10(end(sim) + 1)))

      ## objects to save during simulation --------------------------------------------------------
      times_during <- c(start(sim), end(sim), mod$analysesOutputsTimes) |> unique() |> sort()

      ## fmt: skip
      if (time(sim) %in% times_during) {
        f_rstTimeSinceFire <- file.path(outputPath(sim), paste0("rstTimeSinceFire_year", padYear, ".tif"))
        terra::writeRaster(sim$rstTimeSinceFire, f_rstTimeSinceFire, datatype = "INT2U", overwrite = TRUE)
        sim <- registerOutputs(f_rstTimeSinceFire, sim)

        ## the end is scheduled by init: do not schedule it a second time here
        nextSave <- time(sim) + P(sim)$summaryInterval
        if (time(sim) >= P(sim)$summaryPeriod[1] && time(sim) < P(sim)$summaryPeriod[2] && nextSave != end(sim)) {
          sim <- scheduleEvent(sim, nextSave, "burnSummaries", "save_single", .last())
        }
      }

      ## objects to save at end of simulation -----------------------------------------------------
      if (time(sim) == end(sim)) {
        f_burnMap <- file.path(outputPath(sim), paste0("burnMap_year", padYear, ".tif"))
        terra::writeRaster(sim$burnMap, f_burnMap, datatype = "INT2U", overwrite = TRUE)
        sim <- registerOutputs(f_burnMap, sim)

        f_flammableMap <- file.path(outputPath(sim), paste0("flammableMap_year", padYear, ".tif"))
        terra::writeRaster(sim$flammableMap, f_flammableMap, datatype = "INT2U", overwrite = TRUE)
        sim <- registerOutputs(f_flammableMap, sim)
      }
    },
    noEventWarning(sim)
  )
  return(invisible(sim))
}

## event functions ------------------------------------------------------------------

#' Times at which to summarise: `summaryPeriod` is in simulation times, so these are too
#'
#' @param period The resolved `summaryPeriod`.
#' @param interval The `summaryInterval` parameter.
#'
#' @return A numeric vector of simulation times.
summaryTimes <- function(period, interval) seq(period[1], period[2], by = interval)

#' Resolve `summaryPeriod` against the run's times
#'
#' `NA` means the last 30% of the run: years 700 to 1000 of a 1000-year run, and the last 30% of a
#' shorter one (a run longer than 1000 years still uses 700 to 1000). A period outside the run is
#' an error.
#'
#' @param period The `summaryPeriod` parameter, length 2, `NA` for the default.
#' @param start,end Start and end times of the simulation.
#'
#' @return An integer vector of length 2.
resolveSummaryPeriod <- function(period, start, end) {
  if (anyNA(period)) {
    runLength <- end - start
    period <- if (runLength >= 1000) start + c(700L, 1000L) else c(start + round(0.7 * runLength), end)
    if (runLength < 1000) {
      message("burnSummaries: `summaryPeriod` was not set and the run (", start, " to ", end,
              ") is shorter than 1000 years, so summaries use the last 30% of the run: ",
              period[1], " to ", period[2], ". Set `summaryPeriod` to choose another period.")
    }
  }
  if (min(period) < start || max(period) > end) {
    stop("summaryPeriod ", period[1], " to ", period[2], " is outside the simulation times ",
         start, " to ", end, ". Set summaryPeriod within the run's times.")
  }
  as.integer(period)
}

InitSingle <- function(sim) {
  if (is.null(sim[["flammableMap"]])) {
    ## fireSense uses `flammableRTM`
    stopifnot(suppliedElsewhere("flammableRTM", sim))

    sim <- objectSynonyms(sim, list(c("flammableRTM", "flammableMap")))
  }

  if (is.null(sim$rstTimeSinceFire)) {
    ## fireSense uses `nonForest_timeSinceDisturbance`
    if (!is.null(sim$nonForest_timeSinceDisturbance)) {
      sim$rstTimeSinceFire <- reproducible::postProcess(
        sim$nonForest_timeSinceDisturbance,
        to = sim$flammableMap
      )
    } else {
      sim$rstTimeSinceFire <- LandR::prepInputsStandAgeMap(
        dataSource = "SCANFI",
        dataYear = P(sim)$dataYear,
        ageFun = "terra::rast",
        cropTo = sim$flammableMap,
        maskTo = sim$flammableMap,
        destinationPath = inputPath(sim)
      )

      ## non-flammable areas are permanent
      sim$rstTimeSinceFire[sim$flammableMap[] == 0L] <- NA
      sim$rstTimeSinceFire[] <- as.integer(sim$rstTimeSinceFire[])
    }
  }

  ## sanity check
  terra::compareGeom(
    sim$rstCurrentBurn,
    sim$flammableMap,
    sim$rstTimeSinceFire,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE
  )

  return(invisible(sim))
}

InitMulti <- function(sim) {
  ## check for necessary output files -----------------------------------------------
  mod$useOutputs <- NROW(sim$outputsDF) > 0
  mod$outputFiles <- if (mod$useOutputs) unique(sim$outputsDF$file)
  mod$repDirs <- if (mod$useOutputs) {
    ## the replicate directories as they are named on disk (e.g. rep1 or rep01)
    dirnamesFromSet(mod$outputFiles, P(sim)$reps, leafOnly = FALSE)
  } else {
    file.path(outputPath(sim), sprintf("rep%02d", P(sim)$reps))
  }
  mod$allReps <- basename(mod$repDirs)

  ## assigned back: P(sim)$simTimes is read downstream, not just for padding
  P(sim)$simTimes <- resolveSimYears(P(sim)$simTimes, sim)
  pad <- padYears(P(sim)$simTimes)

  ## all reps have same flammable map
  flm <- repFiles(mod$repDirs, paste0("flammableMap_year", pad$end, ".tif"), mod$outputFiles)[1]

  stopifnot(file.exists(flm))

  flammableMap <- terra::rast(flm)
  pixelSize <- terra::res(flammableMap) ## keep both x and y dimensions

  fbms <- repFiles(mod$repDirs, paste0("burnMap_year", pad$end, ".tif"), mod$outputFiles)
  burnMaps <- Map(rep = mod$allReps, fbm = fbms, function(rep, fbm) {
    message(paste("Loading burn maps for rep", rep, "..."))

    stopifnot(file.exists(fbm))

    cumulBurnMap <- terra::rast(fbm)

    ## sanity check
    terra::compareGeom(cumulBurnMap, flammableMap, res = TRUE)

    ## mean annual cumulative burn map
    cumulBurnMap / (P(sim)$simTimes[2] - P(sim)$simTimes[1])
  }) |>
    terra::rast() |>
    terra::app(sum, na.rm = TRUE)

  meanAnnualCumulBurnMap <- burnMaps / length(mod$allReps)

  ## The cells both maps cover, so that the historic and simulated rates are means over the same
  ## cells: flammable, inside `studyArea`, with a fire return interval (see burnableCells()).
  burnable <- burnableCells(flammableMap, sim$studyArea, sim$fireReturnInterval)

  firePolys <- if (exists("firePolys", envir(sim))) {
    ## polygons supplied on the simList: use them rather than re-fetching
    sim$firePolys |>
      tidyterra::bind_spat_rows() |>
      tidyterra::mutate(YEAR = as.integer(YEAR))
  } else {
    ## Observed fire perimeters: NBAC (National Burned Area Composite -- satellite-derived,
    ## 1972-present, the preferred source) supplemented with NFDB polygons ONLY for the years
    ## before NBAC begins (historicFirePolys()). Older NFDB perimeters are aerial sketches that
    ## overestimate burned area, so NBAC is authoritative wherever it exists. Both are national CWFIS
    ## downloads, harmonised (tolerant YEAR/SIZE_HA columns) + clipped to the sim grid via
    ## fireregimetools::fetch_nbac_polys() / fetch_nfdb_polys().
    message("preparing historical cumulative burn map using NBAC perimeters (+ NFDB backfill)...")
    dst <- inputPath(sim)

    ## fireregimetools fetches + harmonises both archives: each is downloaded once per `dest`
    ## (with R's 60 s timeout raised -- it silently truncates a 1.2 GB transfer -- and staged
    ## through a `.part` file), and an extraction is verified file-by-file against the archive
    ## manifest before it is reused. That guard is not optional here: an interrupted extraction
    ## left the 1.88 GB NBAC .shp at 567 MB, GDAL logged 74,178 read errors, and sf::st_read()
    ## still returned the full feature count (the .shx index was intact), so these historic
    ## summaries were once built from ~30% of the record with nothing failing.
    ## fire_years is left at its default (NULL): keep every year the records contain.
    nbac <- fireregimetools::fetch_nbac_polys(flammableMap, dest = dst)
    nfdb <- fireregimetools::fetch_nfdb_polys(flammableMap, dest = dst)

    historicFirePolys(nbac, nfdb, nbacStartYear = P(sim)$nbacStartYear)
  }

  ## the `historicYears` window, and no undated fire (`YEAR` 0, -9999 or missing; historicFirePolys()
  ## already drops these from the downloaded record, but supplied polygons may hold them)
  hy <- P(sim)$historicYears
  firePolys <- firePolys[which(firePolys$YEAR > 0 & (anyNA(hy) | (firePolys$YEAR >= hy[1] & firePolys$YEAR <= hy[2]))), ]
  recordYears <- historicRecordYears(firePolys, flammableMap, burnable, hy)
  meanAnnualCumulBurnMapHistoric <- historicBurnMap(firePolys, flammableMap, burnable, recordYears)
  ## a burnable cell that no replicate burned has a rate of 0
  meanAnnualCumulBurnMap <- onCells(meanAnnualCumulBurnMap, burnable)

  ## Both maps are written, so the two can be compared cell by cell outside the module.
  maps <- list(
    meanAnnualCumulBurnMapHistoric = meanAnnualCumulBurnMapHistoric,
    meanAnnualCumulBurnMap = meanAnnualCumulBurnMap
  )
  for (nm in names(maps)) {
    f <- file.path(outputPath(sim), paste0(nm, ".tif"))
    terra::writeRaster(maps[[nm]], f, overwrite = TRUE)
    sim <- registerOutputs(f, sim)
  }
  message(sprintf(
    "  mean annual burn rate in the %d burnable cells (share burned per year): historic %s, simulated %s",
    sum(burnable),
    format(mean(terra::values(meanAnnualCumulBurnMapHistoric, mat = FALSE)[burnable]), digits = 4),
    format(mean(terra::values(meanAnnualCumulBurnMap, mat = FALSE)[burnable]), digits = 4)
  ))
  message("...done")

  ## 1 on the burnable cells, NA elsewhere: the land of the observed fire regime, as of both maps
  mod$burnableMap <- terra::rast(flammableMap, vals = ifelse(burnable, 1L, NA_integer_))
  mod$firePolys <- firePolys
  mod$recordYears <- recordYears
  mod$meanAnnualCumulBurnMap <- meanAnnualCumulBurnMap
  mod$meanAnnualCumulBurnMapHistoric <- meanAnnualCumulBurnMapHistoric
  mod$pixelSize <- pixelSize

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

## The file `fname` in each replicate directory. With `files` (the files listed in `outputsDF`),
## the directory and file name must match exactly, so that `rep1` does not also pick up `rep10`'s
## files; without it, the path is built from the directory.
repFiles <- function(repDirs, fname, files = NULL) {
  if (is.null(files)) {
    return(file.path(repDirs, fname))
  }
  vapply(repDirs, function(d) {
    f <- files[dirname(files) == d & basename(files) == fname]
    if (length(f) != 1L) {
      stop("expected one ", fname, " in ", d, " in outputsDF; found ", length(f))
    }
    f
  }, character(1), USE.NAMES = FALSE)
}

## The historic fire record, with only a `YEAR` column: NBAC perimeters from `nbacStartYear` on, and
## NFDB polygons for the years before it. The boundary is a fixed year, not the years NBAC has inside
## the study area: a year from 1972 on with no NBAC fire here is a year in which NBAC mapped none, and
## filling it from the NFDB would mix the two products. NFDB polygons with no year (`YEAR` -9999, 0
## or missing) are dropped.
historicFirePolys <- function(nbac, nfdb, nbacStartYear) {
  nbac <- nbac[which(nbac$YEAR >= nbacStartYear), "YEAR"]
  nDated <- sum(nfdb$YEAR > 0, na.rm = TRUE)
  backfill <- nfdb[which(nfdb$YEAR > 0 & nfdb$YEAR < nbacStartYear), "YEAR"]

  years <- sort(unique(backfill$YEAR))
  message(sprintf(
    "  NBAC from %d: %d polygon(s); backfilling %d NFDB-only year(s): %s",
    nbacStartYear, nrow(nbac), length(years), paste(years, collapse = ", ")
  ))
  if (nrow(backfill) < nrow(nfdb)) {
    message(sprintf(
      "  NFDB polygons not used: %d with no year, %d from %d on",
      nrow(nfdb) - nDated, nDated - nrow(backfill), nbacStartYear
    ))
  }

  if (nrow(backfill) == 0L) {
    nbac
  } else if (nrow(nbac) == 0L) {
    backfill
  } else {
    tidyterra::bind_spat_rows(nbac, backfill)
  }
}

## The cells both burn maps cover (multi mode), as a logical vector over the cells of `flammableMap`:
## flammable (not NA, not 0), with the centre inside `studyArea`, and with a fire return interval
## (not NA, above 0). For LandMine these are the cells where it ignites fires: it masks its spread
## and rate-of-spread maps to `studyArea` and draws ignitions only where there is an interval. Its own
## FRI diagnostics use the same cells. Its fires can also spread into flammable cells of the study
## area that have no interval; those cells are left out here, as they are in the diagnostics.
##
## Through 1.1.0.9002 the historic map covered every flammable cell of the grid, which is the study
## area's bounding box, so it took in fires LandMine could not simulate. 14-33% of the grid's
## flammable cells in LandWeb lay outside the study area.
##
## `studyArea` or `fireReturnInterval` may be NULL: a fire model other than LandMine may have
## neither. Each one missing restricts nothing, and a message says so. The maps then cover cells the
## fire model may not burn, but they are still the same cells for both maps.
burnableCells <- function(flammableMap, studyArea = NULL, fireReturnInterval = NULL) {
  fv <- terra::values(flammableMap, mat = FALSE)
  burnable <- !is.na(fv) & fv != 0
  nFlammable <- sum(burnable)
  kept <- "flammable"

  if (is.null(studyArea)) {
    message("  no `studyArea`: the burn maps include any flammable cells outside the area the fire model burns")
  } else {
    sa <- if (inherits(studyArea, "SpatVector")) studyArea else terra::vect(studyArea)
    ## a polygon in another CRS would rasterize to no cell at all
    if (!terra::same.crs(sa, flammableMap)) {
      sa <- terra::project(sa, terra::crs(flammableMap))
    }
    inside <- !is.na(terra::values(terra::rasterize(sa, flammableMap, field = 1L), mat = FALSE))
    burnable <- burnable & inside
    kept <- c(kept, "inside `studyArea`")
  }

  if (is.null(fireReturnInterval)) {
    message("  no `fireReturnInterval`: the burn maps include any flammable cells without a fire return interval")
  } else {
    terra::compareGeom(fireReturnInterval, flammableMap)
    friv <- terra::values(fireReturnInterval, mat = FALSE)
    burnable <- burnable & !is.na(friv) & friv > 0
    kept <- c(kept, "with a fire return interval")
  }

  if (!any(burnable)) {
    stop(
      "No burnable cell: no cell of the grid is ", paste(kept, collapse = ", "), ". ",
      "Check that `studyArea` and `fireReturnInterval` cover the simulated area."
    )
  }
  message(sprintf(
    "  burn maps cover %d burnable cells (%s) of %d flammable cells",
    sum(burnable), paste(kept, collapse = ", "), nFlammable
  ))
  burnable
}

## The years of the observed record (multi mode): `historicYears` when it is set, else the first to the
## last year with a fire in the burnable cells, inclusive. A year in that span with no fire there is a
## year in which nothing burned, so the historic map is a rate per year, as the simulated one is.
## Through 1.1.0.9002 the divisor was the number of distinct years with a fire anywhere in the grid's
## bounding box, which left fire-free years out and so overstated the rate where they are common.
## `integer(0)` when `historicYears` is unset and no fire touches a burnable cell.
historicRecordYears <- function(firePolys, flammableMap, burnable, historicYears = c(NA, NA)) {
  if (!anyNA(historicYears)) {
    years <- seq(as.integer(historicYears[1]), as.integer(historicYears[2]))
  } else {
    years <- integer(0)
    if (nrow(firePolys) > 0L) {
      yearIn <- function(fun) {
        v <- terra::values(terra::rasterize(firePolys, flammableMap, field = "YEAR", fun = fun), mat = FALSE)
        v[burnable]
      }
      firstByCell <- yearIn("min")
      if (any(!is.na(firstByCell))) {
        years <- seq(as.integer(min(firstByCell, na.rm = TRUE)), as.integer(max(yearIn("max"), na.rm = TRUE)))
      }
    }
  }
  if (length(years)) {
    message(sprintf("  historic record %d-%d (%d years)", min(years), max(years), length(years)))
  } else {
    message("  historic record: no fire in the burnable cells; the historic burn map is 0 there")
  }
  years
}

## The historic burn map (multi mode): recorded fires per cell, divided by the years of record
## (`recordYears`, historicRecordYears()), on the cells in `burnable` and NA elsewhere. A fire counts
## in a cell when its polygon covers the cell's centre, and a cell burned twice counts twice.
historicBurnMap <- function(firePolys, flammableMap, burnable, recordYears) {
  counts <- if (nrow(firePolys) > 0L) {
    terra::rasterize(firePolys, flammableMap, field = "YEAR", fun = "count", background = 0)
  } else {
    terra::rast(flammableMap, vals = 0)
  }
  onCells(if (length(recordYears)) counts / length(recordYears) else counts, burnable)
}

## `x` on the cells in `burnable` only: NA elsewhere, and 0 where `x` is NA inside them (no burn
## recorded is a burn rate of 0).
onCells <- function(x, burnable) {
  v <- terra::values(x, mat = FALSE)
  v[is.na(v)] <- 0
  v[!burnable] <- NA
  terra::setValues(x, v)
}

FireSummaries <- function(sim) {
  pad <- padYears(P(sim)$simTimes)

  studyAreaName <- P(sim)$.studyAreaName

  if (mod$useOutputs) {
    ## registerOutputs() records the CSV; the parquet partition is not registered, so
    ## outputs(sim) is what drives this arm.
    f_csvs <- repFiles(mod$repDirs, "burnSummaries_fireSizes.csv", mod$outputFiles)
    sim$fireSizes <- lapply(f_csvs, function(f_fireSizes) {
      stopifnot(file.exists(f_fireSizes))

      data.table::fread(f_fireSizes)
    }) |>
      data.table::rbindlist()
  } else {
    ## read every replicate's fire-size parquet partition as one lazy Arrow dataset;
    ## open_burn_dataset() skips reps with no fires / missing output rather than erroring.
    roots <- file.path(mod$repDirs, "burnSummaries_fireSizes")
    ds <- fireregimetools::open_burn_dataset(roots)
    sim$fireSizes <- if (is.null(ds)) {
      ## No replicate had a fire, so none wrote a parquet (write_burn_parquet() skips an empty
      ## table) and there is no schema to read. A table with no columns wrote an empty CSV and
      ## stopped plotFun() at `simArea`, so give it the columns and types a replicate with
      ## LandMine fires reads back: create_fireSizes' columns, plus the parquet's `replicate`.
      data.table::data.table(
        simArea = character(0),
        rep = integer(0),
        year = character(0),
        simSize = integer(0),
        expSize = integer(0),
        fireID = integer(0),
        attempt = integer(0),
        targetSize = numeric(0),
        replicate = integer(0)
      )
    } else {
      data.table::as.data.table(dplyr::collect(ds))
    }
  }

  f_out <- file.path(outputPath(sim), paste0("burnSummaries_fireSizes_allReps.csv"))
  data.table::fwrite(sim$fireSizes, f_out)

  ## TODO: add this file to list of outputs
  sim <- registerOutputs(f_out, sim)

  ## the fire regime per replicate and observed, and the fire per block of `summaryInterval` years
  pixelHa <- prod(mod$pixelSize) / 10^4
  repIDs <- as.integer(sub("^\\D*", "", mod$allReps))
  mod$fireRegime <- rbind(
    simulatedFireRegime(sim$fireSizes, repIDs, P(sim)$simTimes, pixelHa),
    observedFireRegime(mod$firePolys, mod$meanAnnualCumulBurnMapHistoric, mod$burnableMap, mod$recordYears)
  )
  mod$fireBlocks <- fireBlockTable(sim$fireSizes, repIDs, P(sim)$simTimes, P(sim)$summaryInterval, pixelHa)
  for (x in c("fireRegime", "fireBlocks")) {
    f <- file.path(outputPath(sim), paste0("burnSummaries_", x, ".csv"))
    data.table::fwrite(mod[[x]], f)
    sim <- registerOutputs(f, sim)
  }

  return(invisible(sim))
}

## The simulated fire regime, one row per replicate in `reps` (a replicate with no fire has 0 fires):
## the number of fires, the area burned per year (all fires' cells over the years from `simTimes[1]` to
## `simTimes[2]`, the same years as the simulated mean annual burn map) and the largest fire, in ha.
## `fireSizes` has columns rep and simSize (cells).
simulatedFireRegime <- function(fireSizes, reps, simTimes, pixelHa) {
  nYears <- simTimes[2] - simTimes[1]
  do.call(rbind, lapply(reps, function(r) {
    s <- as.numeric(fireSizes$simSize[fireSizes$rep == r])
    s <- s[s > 0]
    data.frame(source = "simulated", rep = r, firstYear = simTimes[1], lastYear = simTimes[2],
               nFires = length(s), areaBurnedHaPerYear = sum(s) * pixelHa / nYears,
               largestFireHa = if (length(s)) max(s) * pixelHa else 0)
  }))
}

## The observed fire regime over the same land as the simulated one: area burned per year from the historic
## mean annual burn map (`historic`, NA off the burnable cells), and the largest fire as the burnable land
## inside its polygon (the simulation burns only those cells). `burnableMap` is 1 on the burnable cells
## (burnableCells()), NA elsewhere. `recordYears` is historicRecordYears(); with no record, the first and
## last years are NA.
observedFireRegime <- function(firePolys, historic, burnableMap, recordYears) {
  pixelHa <- prod(terra::res(burnableMap)) / 10^4
  largest <- if (nrow(firePolys)) {
    max(terra::extract(burnableMap, firePolys, fun = sum, na.rm = TRUE, ID = FALSE)[[1]], na.rm = TRUE) * pixelHa
  } else 0
  span <- if (length(recordYears)) range(recordYears) else c(NA_integer_, NA_integer_)
  data.frame(source = "observed", rep = NA_integer_, firstYear = span[1], lastYear = span[2],
             nFires = nrow(firePolys),
             areaBurnedHaPerYear = sum(terra::values(historic, mat = FALSE), na.rm = TRUE) * pixelHa,
             largestFireHa = if (is.finite(largest)) largest else 0)
}

## Fire per block of `blockYears` years per replicate: number of fires and area burned (ha). A block with
## no fire is 0, not missing. `block` = the block's first year; blocks start at `simTimes[1]`, and a last
## block shorter than `blockYears` is left out.
fireBlockTable <- function(fireSizes, reps, simTimes, blockYears, pixelHa) {
  starts <- seq(simTimes[1], simTimes[2] - blockYears, by = blockYears)
  yr <- as.numeric(fireSizes$year)
  sz <- as.numeric(fireSizes$simSize)
  out <- expand.grid(rep = reps, block = starts)
  inBlock <- lapply(seq_len(nrow(out)), function(i)
    fireSizes$rep == out$rep[i] & yr >= out$block[i] & yr < out$block[i] + blockYears & sz > 0)
  out$nFires <- vapply(inBlock, sum, numeric(1))
  out$areaBurnedHa <- vapply(inBlock, function(k) sum(sz[k]) * pixelHa, numeric(1))
  out
}

## Area burned and number of fires per block, one line per replicate, with the observed burn rate scaled
## to a block (dashed) on the area panel. `data` = list(blocks = fireBlockTable(), observed = one row of
## observedFireRegime(), blockYears =, title =).
fireBlockPlot <- function(data) {
  b <- data$blocks
  long <- rbind(
    data.frame(rep = b$rep, block = b$block, value = b$areaBurnedHa / 1000,
               panel = sprintf("Area burned per %d years (thousand ha)", data$blockYears)),
    data.frame(rep = b$rep, block = b$block, value = b$nFires,
               panel = sprintf("Number of fires per %d years", data$blockYears))
  )
  obs <- data.frame(panel = sprintf("Area burned per %d years (thousand ha)", data$blockYears),
                    y = data$observed$areaBurnedHaPerYear * data$blockYears / 1000)
  ggplot2::ggplot(long, ggplot2::aes(block + data$blockYears / 2, value, colour = factor(rep))) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_hline(data = obs, ggplot2::aes(yintercept = y), linetype = "dashed", colour = "black") +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Simulation year (middle of the block)", y = NULL, colour = "Replicate",
                  title = data$title,
                  caption = sprintf("dashed: observed %d-%d burn rate, scaled to %d years (same land, flammable cells)",
                                    data$observed$firstYear, data$observed$lastYear, data$blockYears)) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

## Side-by-side historic and simulated mean annual burn maps on one fill scale.
cumulBurnMapPlots <- function(historic, simulated, studyAreaName, use_palette = "muted") {
  limits <- range(c(terra::values(historic), terra::values(simulated)), finite = TRUE)
  legendTitle <- "Mean annual burn\n(burns per cell per year)"

  mapPlot <- function(r, title) {
    ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = r) +
      tidyterra::scale_fill_whitebox_c(palette = use_palette, limits = limits, name = legendTitle) +
      ggplot2::theme_bw() +
      ggspatial::annotation_north_arrow(
        location = "bl",
        which_north = "true",
        pad_x = ggplot2::unit(0.25, "in"),
        pad_y = ggplot2::unit(0.25, "in"),
        style = north_arrow_fancy_orienteering
      ) +
      ggplot2::xlab("Longitude") +
      ggplot2::ylab("Latitude") +
      ggplot2::ggtitle(paste(title, "mean annual cumulative burn map for", studyAreaName))
  }

  (mapPlot(historic, "Historic") | mapPlot(simulated, "Simulated")) +
    patchwork::plot_layout(guides = "collect")
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  studyAreaName <- P(sim)$.studyAreaName

  ## cumulative burn maps
  use_palette = "muted" # "bl_yl_rd"

  ggCumulBurnMap <- cumulBurnMapPlots(
    mod$meanAnnualCumulBurnMapHistoric,
    mod$meanAnnualCumulBurnMap,
    studyAreaName,
    use_palette
  )

  if ("png" %in% P(sim)$.plots) {
    fggCumulBurnMap <- file.path(figurePath(sim), "cumulative_burn_maps.png")
    ggplot2::ggsave(fggCumulBurnMap, ggCumulBurnMap, height = 10, width = 20, type = "cairo")
    sim <- registerOutputs(fggCumulBurnMap, sim)
  }

  ## fire per block of summaryInterval years, by replicate, against the observed rate
  Plots(
    data = list(blocks = mod$fireBlocks, observed = mod$fireRegime[mod$fireRegime$source == "observed", ],
                blockYears = P(sim)$summaryInterval, title = paste("Fire per", P(sim)$summaryInterval, "years in", studyAreaName)),
    fn = fireBlockPlot, filename = "fire_per_block", path = figurePath(sim),
    types = P(sim)$.plots, ggsaveArgs = list(width = 8, height = 6, units = "in")
  )

  ## fire size histograms w/ median fire sizes
  pixelSizeHa <- prod(mod$pixelSize) / 10^4

  subsetDT <- sim$fireSizes[simArea == studyAreaName & (expSize > 0 | simSize > 0), ]

  ## scfm and fireSense don't set target fire sizes, but LandMine does
  fireModelUsesTargetSize <- isFALSE(all(is.na(subsetDT$expSize)))

  if (isTRUE(fireModelUsesTargetSize)) {
    subsetDT[, expSizeHa := expSize * pixelSizeHa]
  }

  subsetDT[, simSizeHa := simSize * pixelSizeHa]

  ## fire-size distribution (per Dave's email / Eliot): number of fires per log-size bin, with the
  ## median log fire size per bin overlaid on a secondary axis. fireregimetools::fire_size_histogram()
  ## takes raw sizes (ha) and logs internally, replacing the bespoke hist() + stat_summary_bin() code.
  ggHistSim <- fireregimetools::fire_size_histogram(
    subsetDT,
    size_col = "simSizeHa",
    size_unit = "ha",
    title = paste("Total simulated number and size of fires in", studyAreaName)
  )

  if (isTRUE(fireModelUsesTargetSize)) {
    ggHistExp <- fireregimetools::fire_size_histogram(
      subsetDT,
      size_col = "expSizeHa",
      size_unit = "ha",
      title = paste("Total expected number and size of fires in", studyAreaName)
    )
  }

  ## fire_size_histogram() returns NULL when no fire has a positive size (e.g. every replicate
  ## was fire-free), so there is no histogram to save; skip it, as the expected-size plots are
  ## skipped for a fire model without target sizes.
  if ("png" %in% P(sim)$.plots) {
    if (isTRUE(fireModelUsesTargetSize) && !is.null(ggHistExp)) {
      fggHistExp <- file.path(figurePath(sim), "expected_number_size_fires.png")
      ggplot2::ggsave(fggHistExp, ggHistExp, height = 10, width = 10, type = "cairo")
      sim <- registerOutputs(fggHistExp, sim)
    }

    if (is.null(ggHistSim)) {
      message("no simulated fires in ", studyAreaName, "; skipping the fire-size histogram")
    } else {
      fggHistSim <- file.path(figurePath(sim), "simulated_number_size_fires.png")
      ggplot2::ggsave(fggHistSim, ggHistSim, height = 10, width = 10, type = "cairo")
      sim <- registerOutputs(fggHistSim, sim)
    }
  }

  if (isTRUE(fireModelUsesTargetSize)) {
    ## exp vs sim fire sizes
    ggExpVsSim <- ggplot2::ggplot(subsetDT, ggplot2::aes(x = expSizeHa, y = simSizeHa)) +
      ggplot2::geom_smooth(method = lm) +
      ggplot2::scale_x_continuous(limits = c(0, NA)) +
      ggplot2::scale_y_continuous(limits = c(0, NA)) +
      ggplot2::xlab("Expected fire size (ha)") +
      ggplot2::ylab("Simulated fire size (ha)") +
      ggplot2::ggtitle(paste("Expected vs. simulated fire sizes in", studyAreaName)) +
      ggplot2::theme_bw() +
      ggplot2::geom_abline(slope = 1, lty = "dotted")

    ggExpVsSimHex <- ggplot2::ggplot(subsetDT, ggplot2::aes(x = expSizeHa, y = simSizeHa)) +
      ggplot2::geom_hex(bins = 50) +
      ggplot2::xlab("Expected fire size (ha)") +
      ggplot2::ylab("Simulated fire size (ha)") +
      ggplot2::ggtitle(paste("Expected vs. simulated fire sizes in", studyAreaName)) +
      ggplot2::theme_bw() +
      ggplot2::geom_abline(slope = 1, lty = "dotted")

    if ("png" %in% P(sim)$.plots) {
      ## NOTE: keep 1:1 aspect ratio on these plots
      fggExpVsSim <- file.path(figurePath(sim), "exp_vs_sim_fire_sizes.png")
      ggplot2::ggsave(
        filename = fggExpVsSim,
        plot = ggExpVsSim,
        height = 10,
        width = 10,
        type = "cairo"
      )
      sim <- registerOutputs(fggExpVsSim, sim)

      fggExpVsSimHex <- file.path(figurePath(sim), "exp_vs_sim_fire_sizes_hex.png")
      ggplot2::ggsave(
        filename = fggExpVsSimHex,
        plot = ggExpVsSimHex,
        height = 10,
        width = 10,
        type = "cairo"
      )
      sim <- registerOutputs(fggExpVsSimHex, sim)
    }
  }

  ## TODO: is it worth testing fire size distributions? (very slow, and plots show they're bang-on)
  # kSamples::ad.test(subsetDT$simSize, subsetDT$expSize) ## TODO: output this somewhere...

  ## TODO: registerOutputs
  # registerOutputs(c(fggCumulBurnMap, fggHistExp, fggHistSim, fggExpVsSim, fggExpVsSimHex), sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  ## the same data year as the vegetation modules (Biomass_core, Biomass_borealDataPrep, ...)
  paramCheckOtherMods(sim, "dataYear", ifSetButDifferent = "error")
  ## NOTE: fireSize table can be derived from 'burnSummary' or 'fireSize' objects,
  ## so don't implement a hard requirement for either here.
  if (P(sim)$mode == "single") {
    stopifnot(suppliedElsewhere("burnMap", sim))

    ## Create the initial rstTimeSinceFire HERE (in addition to InitSingle) when
    ## flammableMap is supplied directly, so modules that init BEFORE burnSummaries --
    ## e.g. LandMine's Init compareGeom() -- see a non-NULL rstTimeSinceFire (this module
    ## supersedes timeSinceFire, whose .inputObjects used to create it). fireSense supplies
    ## flammableMap via the flammableRTM synonym in InitSingle, so this block is skipped
    ## there (flammableMap not yet present) and InitSingle creates it as before.
    if (is.null(sim$rstTimeSinceFire) && !is.null(sim[["flammableMap"]])) {
      if (!is.null(sim$nonForest_timeSinceDisturbance)) {
        sim$rstTimeSinceFire <- reproducible::postProcess(
          sim$nonForest_timeSinceDisturbance,
          to = sim$flammableMap
        )
      } else {
        sim$rstTimeSinceFire <- LandR::prepInputsStandAgeMap(
          rasterToMatch = sim$flammableMap, ## resample SCANFI (30 m) to the flammableMap grid
          dataSource = "SCANFI",
          dataYear = P(sim)$dataYear,
          ageFun = "terra::rast",
          maskWithRTM = TRUE,
          destinationPath = inputPath(sim)
        )

        ## non-flammable areas are permanent
        sim$rstTimeSinceFire[sim$flammableMap[] == 0L] <- NA
        sim$rstTimeSinceFire[] <- as.integer(sim$rstTimeSinceFire[])
      }
    }
  } else if (P(sim)$mode == "multi") {
    stopifnot(suppliedElsewhere("reportingPolygons", sim))
  }

  return(invisible(sim))
}
