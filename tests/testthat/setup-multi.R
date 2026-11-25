## Fixtures for multi-mode tests: a stand-in for LandMine, single-mode replicates that write what multi mode
## reads, and a multi-mode run over them with the national fire records mocked.
##
## A setup file, not a helper: these functions read `testPaths` (setup.R) when called. pkgload sources
## helper-*.R files into the package namespace, which cannot see objects setup files create in the test
## environment, so as helpers they stopped with "object 'testPaths' not found".

## A stand-in for LandMine: one `fireSizes` table per burn year with LandMine's columns (fireID,
## attempt, targetSize, size, maxSize), zero rows in a fire-free year. With `fires = TRUE` it burns
## cell 1 every year; with `FALSE`, nothing.
writeLandMineStandIn <- function(dir) {
  d <- file.path(dir, "landMineStandIn")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "defineModule(sim, list(",
      '  name = "landMineStandIn", description = "test stand-in for LandMine fire outputs",',
      '  keywords = "test", authors = person("Test", "Author"), childModules = character(0),',
      '  version = list(landMineStandIn = "0.0.1"), timeframe = as.POSIXlt(c(NA, NA)),',
      '  timeunit = "year", citation = list(), documentation = list(),',
      '  reqdPkgs = list("data.table", "terra"),',
      "  parameters = bindrows(",
      '    defineParameter("fires", "logical", TRUE, NA, NA, "burn cell 1 every year")',
      "  ),",
      "  inputObjects = bindrows(",
      '    expectsInput("flammableMap", "SpatRaster", "flammable cells")',
      "  ),",
      "  outputObjects = bindrows(",
      '    createsOutput("burnMap", "SpatRaster", "cumulative burns"),',
      '    createsOutput("fireSizes", "list", "fire sizes by year"),',
      '    createsOutput("rstCurrentBurn", "SpatRaster", "cells burned this year")',
      "  )",
      "))",
      "doEvent.landMineStandIn <- function(sim, eventTime, eventType) {",
      "  switch(eventType,",
      "    init = {",
      "      sim$rstCurrentBurn <- terra::rast(sim$flammableMap, vals = 0L)",
      "      sim$burnMap <- terra::rast(sim$flammableMap, vals = 0L)",
      "      sim$fireSizes <- list()",
      '      sim <- scheduleEvent(sim, start(sim) + 1, "landMineStandIn", "burn", 2.5)',
      "    },",
      "    burn = {",
      "      n <- if (isTRUE(P(sim)$fires)) 1L else 0L",
      "      sim$rstCurrentBurn[] <- 0L",
      "      if (n > 0L) sim$rstCurrentBurn[1] <- 1L",
      "      sim$burnMap <- sim$burnMap + sim$rstCurrentBurn",
      "      sim$fireSizes[[as.character(time(sim))]] <- data.table::data.table(",
      "        fireID = seq_len(n), attempt = rep(1L, n), targetSize = rep(1, n),",
      "        size = rep(1L, n), maxSize = rep(1L, n)",
      "      )",
      '      sim <- scheduleEvent(sim, time(sim) + 1, "landMineStandIn", "burn", 2.5)',
      "    }",
      "  )",
      "  invisible(sim)",
      "}"
    ),
    file.path(d, "landMineStandIn.R")
  )
  dir
}

testGrid <- function() {
  terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = "EPSG:3400")
}

## One replicate in single mode: writes what multi mode reads (fire-size CSV and parquet, and the
## end-of-run burnMap and flammableMap).
runSingleRep <- function(out, rep, fires, mods, end = 4) {
  r <- testGrid()
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = end),
    params = list(
      landMineStandIn = list(fires = fires),
      burnSummaries = list(
        mode = "single", reps = rep, summaryPeriod = c(2L, 4L), summaryInterval = 1L,
        .studyAreaName = "test"
      )
    ),
    modules = c("landMineStandIn", "burnSummaries"),
    loadOrder = c("landMineStandIn", "burnSummaries"),
    objects = list(
      flammableMap = terra::rast(r, vals = c(1L, 1L, 1L, 0L)),
      rstTimeSinceFire = terra::rast(r, vals = c(10L, 20L, 30L, NA))
    ),
    paths = utils::modifyList(
      testPaths,
      list(modulePath = c(testPaths$modulePath, mods), outputPath = out)
    )
  )
}

## Multi mode over the replicates under `out`. The national fire records are replaced by one
## polygon over the grid, so nothing is downloaded. With `outputsDF` (a table with a `file` column) multi
## mode reads the files it lists, as LandWeb runs it, rather than the replicate directories under `out`.
runMulti <- function(out, reps, end = 4, historicYears = c(NA_integer_, NA_integer_),
                     regimePeriod = c(NA_integer_, NA_integer_), outputsDF = NULL, plots = "png",
                     objects = list()) {
  polys <- terra::as.polygons(terra::ext(testGrid()), crs = "EPSG:3400")
  polys$YEAR <- 2000L
  polys$FIRECAUS <- "Natural" # the loaders keep each source's cause column (NBAC FIRECAUS, NFDB CAUSE)
  polys$CAUSE <- "N"
  testthat::local_mocked_bindings(
    fetch_nbac_polys = function(...) polys,
    fetch_nfdb_polys = function(...) polys,
    .package = "fireregimetools"
  )
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = end),
    params = list(burnSummaries = list(
      mode = "multi", reps = reps, simTimes = c(0, end), summaryPeriod = c(2L, 4L),
      summaryInterval = 1L, minCellHa = 1, .studyAreaName = "test", .plots = plots, historicYears = historicYears,
      regimePeriod = regimePeriod
    )),
    modules = "burnSummaries",
    objects = utils::modifyList(objects, list(outputsDF = outputsDF)), ## a NULL outputsDF adds nothing
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )
}

## The files multi mode reads from each replicate directory under `out`, as `outputsDF` lists them.
repOutputs <- function(out, reps, end = 4) {
  dirs <- file.path(out, sprintf("rep%02d", reps))
  data.table::data.table(file = normalizePath(c(
    file.path(dirs, sprintf("flammableMap_year%d.tif", end)),
    file.path(dirs, sprintf("burnMap_year%d.tif", end)),
    file.path(dirs, "burnSummaries_fireSizes.csv")
  ), winslash = "/", mustWork = TRUE))
}
