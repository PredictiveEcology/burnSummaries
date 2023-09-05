defineModule(sim, list(
  name = "burnSummaries",
  description = "",
  keywords = "",
  authors = c(
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre"),
           comment = c(ORCID = "0000-0001-7146-8135"))
  ),
  childModules = character(0),
  version = list(burnSummaries = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "burnSummaries.Rmd"),
  reqdPkgs = list("data.table", "ggplot2", "kSamples", "LandWebUtils",
                  "patchwork", "raster", "rasterVis", "SpaDES.core"),
  parameters = bindrows(
    defineParameter("reps", "integer", 1L:10L, 1L, NA_integer_,
                    paste("number of replicates/runs per study area.")),
    defineParameter("simOutputPath", "character", outputPath(sim), NA, NA,
                    "Directory specifying the location of the simulation outputs."),
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
    expectsInput("speciesLayers", "RasterStack",
                 desc = "initial percent cover raster layers used for simulation."),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "initial percent cover raster layers used for simulation.")
  ),
  outputObjects = bindrows(
    createsOutput("fireSizes", "data.table", "summary fire sizes table"),
    createsOutput("friSummary", "data.table", "summary fire return interval table")
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
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "burnSummaries", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "burnSummaries", "areaBurned")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "burnSummaries", "fireReturnIntervals")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "burnSummaries", "fireSizes")

      if (isTRUE(P(sim)$upload)) {
        sim <- scheduleEvent(sim, end(sim), "LandWeb_summary", "upload", .last())
      }
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      plotFun(sim) # example of a plotting function
      # schedule future event(s)

      # e.g.,
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "burnSummaries", "plot")

      # ! ----- STOP EDITING ----- ! #
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
    areaBurned = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- AreaBurned(sim)

      # ! ----- STOP EDITING ----- ! #
    },
    fireReturnIntervals = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- FireReturnIntervals(sim)

      # ! ----- STOP EDITING ----- ! #
    },
    fireSizes = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- FireSizes(sim)

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

  flammableMap <- NULL
  lthfc <- NULL
  pixelRes <- NULL

  burnMaps <- lapply(P(sim)$reps, function(rep) {
    fsim <- findSimFile(outputPath(sim), rep)

    tmpSim <- loadSimList(fsim)

    if (rep == 1L) {
      ## all reps have same flammable + LTHFC map
      flammableMap <<- tmpSim[["rstFlammable"]] ## RasterLayer
      lthfc <<- tmpSim[["fireReturnInterval"]] ## RasterLayer
      pixelRes <<- res(tmpSim[["rasterToMatch"]]) ## c(250, 250)
    }

    ## sanity check
    compareRaster(tmpSim[["rstCurrentBurnCumulative"]], tmpSim[["rstFlammable"]], res = TRUE, orig = TRUE)

    ## mean annual cumulative burn map
    tmpSim[["rstCurrentBurnCumulative"]] / (end(tmpSim) - start(tmpSim))
  }) |> raster::stack() |>
    raster::calc(sum, na.rm = TRUE)

  mod$flammmableMap <- flammableMap
  mod$lthfc <- lthfc
  mod$meanAnnualCumulBurnMap <- burnMaps / length(reps)
  mod$pixelRes <- pixelRes

  ## remove non-flammable pixels
  toRm <- which(is.na(flammableMap[]) | flammableMap[] == 0)
  if (length(toRm) > 0) {
    mod$meanAnnualCumulBurnMap[toRm] <- NA
    mod$flammableMap[toRm] <- NA
    mod$lthfc[toRm] <- NA
  }

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

## TODO: this event is LandMine-specific; move there or use conditionally based on existence of fireReturnInterval raster
FireReturnIntervals <- function(sim) {
  friValues <- getValues(mod$lthfc) |>
    unique() |>
    na.omit() |>
    sort()

  simFRIs <- vapply(friValues, function(fri) {
    pixIds <- which(getValues(mod$lthfc) == fri)
    1 / (sum(mod$meanAnnualCumulBurnMap[pixIds]) / (length(pixIds)))
  }, numeric(1))

  sim$friSummary <- data.table(
    simArea = P(sim)$.studyAreaName,
    LTHFC = friValues,
    FRI = simFRIs,
    stringsAsFactors = FALSE
  )

  f <- file.path(outputPath(sim), paste0("burnSummaries_FRI_table.csv"))
  fwrite(sim$friSummary, f) ## TODO: add this file to list of outputs

  return(invisible(sim))
}

FireSizes <- function(sim) {
  sim$fireSizes <- lapply(P(sim)$reps, function(rep) {
    fsim <- findSimFile(outputPath(sim), rep)

    tmpSim <- loadSimList(fsim)

    if (!is.null(tmpSim[["fireSizes"]])) {
      fs <- rbindlist(tmpSim[["fireSizes"]], idcol = "year")
      fs[, `:=`(simArea = area, rep = rep)]
      setcolorder(fs, c("simArea", "rep", "year", "size", "maxSize"))
      setnames(fs, old = c("size", "maxSize"), new = c("simSize", "expSize"))
    } else {
      NULL
    }
  }) |>
    rbindlist()

  sim$fireSizes[, simArea := P(sim)$.studyAreaName]

  f <- file.path(P(sim)$simOutputPath, paste0("burnSummaries_fireSizes.csv"))
  fwrite(sim$fireSizes, f) ## TODO: add this file to list of outputs

  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
browser()
  ## TODO: use Plots

  pixelSizeHa <- prod(mod$pixelRes) / 10^4

  ## cumulative burn maps
  ggCumulBurnMapExp <- rasterVis::levelplot(
    mod$meanAnnualCumulBurnMapHistoric,
    main = paste("Historic mean annual cumulative burn map for", P(sim)$.studyAreaName),
    margin = FALSE,
    par.settings = magmaTheme ## TODO: use decent colour scheme
  )

  ggCumulBurnMapSim <- rasterVis::levelplot(
    mod$meanAnnualCumulBurnMap,
    main = paste("Simulated mean annual cumulative burn map for", P(sim)$.studyAreaName),
    margin = FALSE,
    par.settings = magmaTheme ## TODO: use decent colour scheme
  )

  ## TODO: plot cumulative burn map from NFDB data beside the simulated cumul burn map
  ggCumulBurnMap <- ggCumulBurnMapExp + ggCumulBurnMapSim

  fggCumulBurnMap <- file.path(figurePath(sim), "cumulative_burn_maps.png")
  png(fggCumulBurnMap, height = 1000, width = 2000)
  print(ggCumulBurnMap)
  dev.off()

  ## LTHFC/FRI polygons
  ggFriPolys <- rasterVis::levelplot(
    mod$lthfc,
    main = paste("Long-term historic fire cycle (LTHFC) map for", P(sim)$.studyAreaName),
    margin = FALSE,
    par.settings = PuOrTheme
  )

  fggFriPolys <- file.path(figurePath(sim), "LandMine_LTHFC_map.png")
  png(fggFriPolys, height = 1000, width = 1000)
  print(ggFriPolys)
  dev.off()

  ## expected vs simulated fire return intervals
  ggFriExpVsSim <- ggplot(friTable, aes(x = LTHFC, y = FRI)) +
    geom_point() +
    xlab("Expected fire return interval (years)") +
    ylab("Simulated fire return interval (years)") +
    ggtitle(paste("Expected vs. simulated fire return intervals in", P(sim)$.studyAreaName)) +
    theme_bw() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    geom_abline(slope = 1, lty = "dotted")

  fggFriExpVsSim <- file.path(figurePath(sim), "burnSummaries_exp_vs_sim_FRI.png")
  ggsave(filename = fggFriExpVsSim,
         plot = ggFriExpVsSim,
         height = 1000, width = 1000) ## NOTE: keep square aspect ratio

  ## fire size histograms w/ medion fire sizes
  subsetDT <- sim$fireSizes[simArea == P(sim)$.studyAreaName & (expSize > 0 | simSize > 0), ]

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

  breaks <- seq(1.0, ceiling(max(maxLogExpSizeHa, maxLogSimSizeHa) / 0.5) * 0.5, 0.5)

  hexp <- hist(subsetDT$logExpSizeHa, breaks = breaks, plot = FALSE)
  hsim <- hist(subsetDT$logSimSizeHa, breaks = breaks, plot = FALSE)

  countsExp <- hexp$counts ## correct counts?
  countsSim <- hsim$counts ## correct counts?

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
    ggtitle(paste("Total expected number and size of fires in", P(sim)$.studyAreaName)) +
    theme_bw() +
    theme(
      axis.title.y.left = element_text(color = y1col),
      axis.text.y.left = element_text(color = y1col),
      axis.title.y.right = element_text(color = y2col),
      axis.text.y.right = element_text(color = y2col)
    )

  fggHistExp <- file.path(figurePath(sim), "simulated_number_size_fires.png")
  ggsave(ggHistExp, fggHistExp, height = 1000, width = 1000)

  ggHistSim <- ggplot(subsetDT, aes(x = logSimSizeHa)) +
    geom_histogram(breaks = breaks, alpha = 0.5, fill = y1col) +
    stat_summary_bin(data = summarySimDT,
                     mapping = aes(x = midsSim, y = medLogSimSizeHa * scaleFactorSim),
                     fun = "identity", geom = "point", breaks = breaks, col = y2col) +
    scale_y_continuous(y1lab, sec.axis = sec_axis(~ . / scaleFactorSim , name = y2lab)) +
    xlab(x_lab) +
    ggtitle(paste("Total simulated number and size of fires in", P(sim)$.studyAreaName)) +
    theme_bw() +
    theme(
      axis.title.y.left = element_text(color = y1col),
      axis.text.y.left = element_text(color = y1col),
      axis.title.y.right = element_text(color = y2col),
      axis.text.y.right = element_text(color = y2col)
    )

  fggHistSim <- file.path(figurePath(sim), "expected_number_size_fires.png")
  ggsave(ggHistSim, fggHistSim, height = 1000, width = 1000)

  ## exp vs sim fire sizes
  ggExpVsSimHex <- ggplot(subsetDT, aes(x = expSizeHa, y = simSizeHa)) +
    geom_hex(bins = 50) +
    xlab("Expected fire size (ha)") +
    ylab("Simulated fire size (ha)") +
    ggtitle(paste("Expected vs. simulated fire sizes in", P(sim)$.sutdyAreaName)) +
    theme_bw() +
    geom_abline(slope = 1, lty = "dotted")

  fggExpVsSimHex <- file.path(figurePath(sim), "burnSummaries_exp_vs_sim_hex.png")
  ggsave(filename = fggExpVsSimHex,
         plot = ggExpVsSimHex,
         height = 1000, width = 1000) ## NOTE: keep square aspect ratio


  ## test fire size distributions (very slow...)
  kSamples::ad.test(subsetDT$simSize, subsetDT$expSize) ## TODO: output this somewhere...

  ##  track which plot files to upload
  mod$files2upload <- c(fggFriPolys, fggFriExpVsSim, fggHistExp, fggHistSim, fggExpVsSimHex)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

## older version of SpaDES.core used here doesn't have this function
if (packageVersion("SpaDES.core" < "2.0.2.9001")) {
  figurePath <- function(sim) {
    file.path(outputPath(sim), "figures", current(sim)[["moduleName"]])
  }
}
