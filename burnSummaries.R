defineModule(sim, list(
  name = "burnSummaries",
  description = "",
  keywords = "",
  authors = c(
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre"),
           comment = c(ORCID = "0000-0001-7146-8135"))
  ),
  childModules = character(0),
  version = list(burnSummaries = "1.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "burnSummaries.Rmd"),
  reqdPkgs = list("data.table", "ggplot2", "kSamples", "patchwork",
                  "rasterVis", "reproducible", "SpaDES.core", "terra"),
  parameters = bindrows(
    defineParameter("reps", "integer", 1L:10L, 1L, NA_integer_,
                    paste("number of replicates/runs per study area.")),
    defineParameter("simOutPrefix", "character", "mySimOut", NA_character_, NA_character_,
                    "saved simList file prefix"),
    defineParameter("simOutputPath", "character", outputPath(sim), NA, NA,
                    "Directory specifying the location of the simulation outputs."),
    defineParameter("simTimes", "numeric", c(NA, NA), NA, NA,
                    "Simulation start and end times when running in 'multi' mode."),
    defineParameter("upload", "logical", FALSE, NA, NA,
                    "if TRUE, uses the `googledrive` package to upload figures."),
    defineParameter("uploadTo", "character", NA, NA, NA,
                    paste("if `upload = TRUE`, a Google Drive folder id corresponding to `.studyAreaName`.")),
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
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("speciesLayers", "SpatRaster",
                 desc = "initial percent cover raster layers used for simulation.")
  ),
  outputObjects = bindrows(
    createsOutput("fireSizes", "data.table", "summary fire sizes table")
  )
))

## event types
#   - type `init` is required for initialization

doEvent.burnSummaries = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim), "burnSummaries", "summary")
      sim <- scheduleEvent(sim, start(sim), "burnSummaries", "plot")

      if (isTRUE(P(sim)$upload)) {
        sim <- scheduleEvent(sim, end(sim), "burnSummaries", "upload", .last())
      }
    },
    summary = {
      sim <- FireSummaries(sim)
    },
    plot = {
      plotFun(sim)
    },
    upload = {
      # ! ----- EDIT BELOW ----- ! #
      browser() ## TODO
      mod$files2upload <- set_names(mod$files2upload, basename(mod$files2upload))

      gid <- as_id(sim$uploadTo[[P(sim)$.studyAreaName]])
      prevUploaded <- drive_ls(gid)
      toUpload <- mod$files2upload[!(basename(mod$files2upload) %in% prevUploaded$name)]
      uploaded <- map(toUpload, ~ drive_upload(.x, path = gid))
      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #
  allReps <- sprintf("rep%02d", P(sim)$reps)
  mod$files2upload <- list()
  flammableMap <- NULL
  pixelRes <- NULL

  burnMaps <- lapply(allReps, function(rep) {
    fsim <- file.path(outputPath(sim), rep,
                      sprintf("%s_%04d.qs", P(sim)$simOutPrefix, P(sim)$simTimes[2]))
    if (!file.exists(fsim)) {
      fsim <- paste0(tools::file_path_sans_ext(fsim), ".rds") ## fallback to rds if qs not used
    }

    stopifnot(file.exists(fsim))

    tmpSim <- suppressMessages(loadSimList(fsim))

    if (rep == allReps[1]) {
      ## all reps have same flammable map
      flammableMap <<- tmpSim[["rstFlammable"]]
      pixelRes <<- res(tmpSim[["rasterToMatch"]])
    }

    ## NOTE: fire models differ in their naming scheme for cumul burn maps
    ## - LandMine uses `rstCurrentBurnCumulative`;
    ## - fireSense uses `burnMap`;
    ## - scfm use `burnMap`;
    cumulBurnMap <- if (!is.null(tmpSim[["rstCurrentBurnCumulative"]])) {
      tmpSim[["rstCurrentBurnCumulative"]]
    } else if (!is.null(tmpSim[["burnMap"]])) {
      tmpSim[["burnMap"]]
    } else {
      stop("not sure which object is used for cumulative burn map")
    }

    ## sanity check
    terra::compareGeom(cumulBurnMap, tmpSim[["rstFlammable"]], res = TRUE)

    ## mean annual cumulative burn map
    cumulBurnMap / (end(tmpSim) - start(tmpSim))
  }) |>
    terra::rast() |>
    terra::app(sum, na.rm = TRUE)

  meanAnnualCumulBurnMap <- burnMaps / length(allReps)

  firePolys <- Cache(
    prepInputs,
    url = "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip",
    fun = "sf::st_read",
    destinationPath = inputPath(sim),
    targetFile = "NFDB_poly_20210707.shp" ## TODO: this shouldn't be needed; needs to be updated
  ) |>
    sf::st_cast("MULTIPOLYGON") |>
    sf::st_transform(crs(flammableMap))

  fireYears <- firePolys[firePolys$YEAR > 0, ][["YEAR"]] |> unique() |> sort()
  meanAnnualCumulBurnMapHistoric <- terra::rasterize(firePolys, flammableMap, field = "YEAR", fun = "count")
  meanAnnualCumulBurnMapHistoric <- meanAnnualCumulBurnMapHistoric / length(fireYears)

  nonFlammable <- which(is.na(terra::values(flammableMap, mat = FALSE)) |
                          terra::values(flammableMap, mat = FALSE) == 0)
  if (length(nonFlammable) > 0) {
    meanAnnualCumulBurnMap[nonFlammable] <- NA
    meanAnnualCumulBurnMapHistoric[nonFlammable] <- NA
    flammableMap[nonFlammable] <- NA
  }

  mod$flammableMap <- flammableMap
  mod$meanAnnualCumulBurnMap <- meanAnnualCumulBurnMap
  mod$meanAnnualCumulBurnMapHistoric <- meanAnnualCumulBurnMapHistoric
  mod$pixelRes <- pixelRes

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

FireSummaries <- function(sim) {
  allReps <- sprintf("rep%02d", P(sim)$reps)

  studyAreaName <- P(sim)$.studyAreaName
  outputDir <- outputPath(sim)

  sim$fireSizes <- lapply(allReps, function(rep) {
    fsim <- file.path(outputPath(sim), rep,
                      sprintf("%s_%04d.qs", P(sim)$simOutPrefix, P(sim)$simTimes[2]))
    if (!file.exists(fsim)) {
      fsim <- paste0(tools::file_path_sans_ext(fsim), ".rds") ## fallback to rds if qs not used
    }

    stopifnot(file.exists(fsim))

    tmpSim <- suppressMessages(loadSimList(fsim))

    if (!is.null(tmpSim[["burnSummary"]])) {
      fs <- copy(tmpSim[["burnSummary"]])
      fs[, `:=`(simArea = studyAreaName,
                size = N, ## use number of pixels (N) instead of areaBurned
                maxSize = NA_integer_, ## NOTE: no "target" or "expected" fire size w/ scfm nor fS
                rep = as.integer(substr(rep, 4, 5)))]
      set(fs, NULL, c("areaBurned", "igLoc", "grp", "PolyID"), NULL)
      setcolorder(fs, c("simArea", "rep", "year", "size", "maxSize"))
      setnames(fs, old = c("size", "maxSize"), new = c("simSize", "expSize"))
    } else if (!is.null(tmpSim[["fireSizes"]])) {
      fs <- rbindlist(tmpSim[["fireSizes"]], idcol = "year")
      fs[, `:=`(simArea = studyAreaName, rep = rep)]
      setcolorder(fs, c("simArea", "rep", "year", "size", "maxSize"))
      setnames(fs, old = c("size", "maxSize"), new = c("simSize", "expSize"))
    } else {
      NULL
    }
  }) |>
    rbindlist()

  f <- file.path(outputDir, paste0("burnSummaries_fireSizes.csv"))
  fwrite(sim$fireSizes, f) ## TODO: add this file to list of outputs

  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  ## TODO: use Plots
  studyAreaName <- P(sim)$.studyAreaName
  pixelSizeHa <- prod(mod$pixelRes) / 10^4

  ## cumulative burn maps
  ggCumulBurnMapExp <- rasterVis::levelplot(
    mod$meanAnnualCumulBurnMapHistoric,
    main = paste("Historic mean annual cumulative burn map for", studyAreaName),
    margin = FALSE,
    par.settings = magmaTheme ## TODO: use decent colour scheme
  )

  ggCumulBurnMapSim <- rasterVis::levelplot(
    mod$meanAnnualCumulBurnMap,
    main = paste("Simulated mean annual cumulative burn map for", studyAreaName),
    margin = FALSE,
    par.settings = magmaTheme ## TODO: use decent colour scheme
  )

  if ("png" %in% P(sim)$.plots) {
    fggCumulBurnMap <- file.path(figurePath(sim), "cumulative_burn_maps.png")
    png(fggCumulBurnMap, height = 1000, width = 2000)
    gridExtra::grid.arrange(ggCumulBurnMapExp, ggCumulBurnMapSim, ncol = 2, nrow = 1)
    dev.off()
  }

  ## fire size histograms w/ median fire sizes
  subsetDT <- sim$fireSizes[simArea == studyAreaName & (expSize > 0 | simSize > 0), ]

  subsetDT[, expSizeHa := expSize * pixelSizeHa]
  subsetDT[, simSizeHa := simSize * pixelSizeHa]

  subsetDT[, logExpSize := log(expSize)]
  subsetDT[, logSimSize := log(simSize)]

  subsetDT[, logExpSizeHa := log(expSize * pixelSizeHa)]
  subsetDT[, logSimSizeHa := log(simSize * pixelSizeHa)]

  ## per Dave's original email:
  ## > What I would like is both the number of disturbances on the y axis,
  ## > and the area of disturbances on a second y-axis graph.
  ## Per Eliot: x-axis uses same bins as histogram, with y-axis of median area burned per bin

  maxLogExpSizeHa <- max(subsetDT$logExpSizeHa)
  maxLogSimSizeHa <- max(subsetDT$logSimSizeHa)

  breaks <- seq(0.0, ceiling(max(maxLogExpSizeHa, maxLogSimSizeHa, na.rm = TRUE) / 0.5) * 0.5, 0.5)

  hexp <- hist(subsetDT$logExpSizeHa, breaks = breaks, plot = FALSE)
  hsim <- hist(subsetDT$logSimSizeHa, breaks = breaks, plot = FALSE)

  countsExp <- hexp$counts
  countsSim <- hsim$counts

  subsetDT[, binIDexp := cut(logExpSizeHa, hexp$breaks)]
  subsetDT[, binIDsim := cut(logSimSizeHa, hsim$breaks)]

  summaryExpDT <- subsetDT[, lapply(.SD, stats::median, na.rm = TRUE), by = binIDexp, .SDcols = "expSizeHa"]
  summarySimDT <- subsetDT[, lapply(.SD, stats::median, na.rm = TRUE), by = binIDsim, .SDcols = "simSizeHa"]
  setnames(summaryExpDT, "expSizeHa", "medExpSizeHa")
  setnames(summarySimDT, "simSizeHa", "medSimSizeHa")

  summaryExpDT <- summaryExpDT[, medLogExpSizeHa := log(medExpSizeHa)]
  summarySimDT <- summarySimDT[, medLogSimSizeHa := log(medSimSizeHa)]

  midsExp <- cbind(
    as.numeric( sub("\\((.+),.*", "\\1", summaryExpDT$binIDexp) ),
    as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", summaryExpDT$binIDexp) )
  ) |>
    rowMeans()
  summaryExpDT <- summaryExpDT[, midsExp := midsExp]

  midsSim <- cbind(
    as.numeric( sub("\\((.+),.*", "\\1", summarySimDT$binIDsim) ),
    as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", summarySimDT$binIDsim) )
  ) |>
    rowMeans()
  summarySimDT <- summarySimDT[, midsSim := midsSim]

  scaleFactorExp <- max(countsExp) / maxLogExpSizeHa
  scaleFactorSim <- max(countsSim) / maxLogSimSizeHa

  y1col <- "grey20"
  y2col <- "darkred"
  y1lab <- "Total number of fires across all simulated years"
  y2lab <- "Median log[fireSize] (ha)"
  x_lab <- "log[fireSize] (ha)"

  ggHistExp <- ggplot(subsetDT, aes(x = logExpSizeHa)) +
    geom_histogram(breaks = breaks, alpha = 0.5, fill = y1col) +
    stat_summary_bin(data = summaryExpDT,
                     mapping = aes(x = midsExp, y = medLogExpSizeHa * scaleFactorExp),
                     fun = "identity", geom = "point", breaks = breaks, col = y2col) +
    scale_y_continuous(y1lab, sec.axis = sec_axis(~ . / scaleFactorExp , name = y2lab)) +
    xlab(x_lab) +
    ggtitle(paste("Total expected number and size of fires in", studyAreaName)) +
    theme_bw() +
    theme(
      axis.title.y.left = element_text(color = y1col),
      axis.text.y.left = element_text(color = y1col),
      axis.title.y.right = element_text(color = y2col),
      axis.text.y.right = element_text(color = y2col)
    )

  ggHistSim <- ggplot(subsetDT, aes(x = logSimSizeHa)) +
    geom_histogram(breaks = breaks, alpha = 0.5, fill = y1col) +
    stat_summary_bin(data = summarySimDT,
                     mapping = aes(x = midsSim, y = medLogSimSizeHa * scaleFactorSim),
                     fun = "identity", geom = "point", breaks = breaks, col = y2col) +
    scale_y_continuous(y1lab, sec.axis = sec_axis(~ . / scaleFactorSim , name = y2lab)) +
    xlab(x_lab) +
    ggtitle(paste("Total simulated number and size of fires in", studyAreaName)) +
    theme_bw() +
    theme(
      axis.title.y.left = element_text(color = y1col),
      axis.text.y.left = element_text(color = y1col),
      axis.title.y.right = element_text(color = y2col),
      axis.text.y.right = element_text(color = y2col)
    )

  if ("png" %in% P(sim)$.plots) {
    fggHistExp <- file.path(figurePath(sim), "expected_number_size_fires.png")
    ggsave(fggHistExp, ggHistExp, height = 10, width = 10, type = "cairo")

    fggHistSim <- file.path(figurePath(sim), "simulated_number_size_fires.png")
    ggsave(fggHistSim, ggHistSim, height = 10, width = 10, type = "cairo")
  }

  ## exp vs sim fire sizes
  ggExpVsSim <- ggplot(subsetDT, aes(x = expSizeHa, y = simSizeHa)) +
    geom_smooth(method = lm) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    xlab("Expected fire size (ha)") +
    ylab("Simulated fire size (ha)") +
    ggtitle(paste("Expected vs. simulated fire sizes in", studyAreaName)) +
    theme_bw() +
    geom_abline(slope = 1, lty = "dotted")

  ggExpVsSimHex <- ggplot(subsetDT, aes(x = expSizeHa, y = simSizeHa)) +
    geom_hex(bins = 50) +
    xlab("Expected fire size (ha)") +
    ylab("Simulated fire size (ha)") +
    ggtitle(paste("Expected vs. simulated fire sizes in", studyAreaName)) +
    theme_bw() +
    geom_abline(slope = 1, lty = "dotted")

  if ("png" %in% P(sim)$.plots) {
    ## NOTE: keep 1:1 aspect ratio on these plots
    fggExpVsSim <- file.path(figurePath(sim), "exp_vs_sim_fire_sizes.png")
    ggsave(filename = fggExpVsSim, plot = ggExpVsSim, height = 10, width = 10, type = "cairo")

    fggExpVsSimHex <- file.path(figurePath(sim), "exp_vs_sim_fire_sizes_hex.png")
    ggsave(filename = fggExpVsSimHex, plot = ggExpVsSimHex, height = 10, width = 10, type = "cairo")
  }

  ## TODO: is it worth testing fire size distributions? (very slow, and plots show they're bang-on)
  # kSamples::ad.test(subsetDT$simSize, subsetDT$expSize) ## TODO: output this somewhere...

  ##  track which plot files to upload
  mod$files2upload <- append(mod$files2upload, list(
    fggCumulBurnMap, fggHistExp, fggHistSim, fggExpVsSim, fggExpVsSimHex
  ))

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

## older version of SpaDES.core used here doesn't have this function
if (packageVersion("SpaDES.core") < "2.0.2.9001") {
  figurePath <- function(sim) {
    file.path(outputPath(sim), "figures", current(sim)[["moduleName"]]) |>
      checkPath(create = TRUE)
  }
}
