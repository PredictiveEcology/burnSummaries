## Multi mode when no replicate had a fire. A fire-free replicate writes no parquet
## (fireregimetools::write_burn_parquet() skips an empty table), so multi mode has no schema to read.
## Through 1.1.0.9000 it built a table with no columns: burnSummaries_fireSizes_allReps.csv was
## written empty ("Input has no columns") and plotFun() stopped with "Object 'simArea' not found".

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
## polygon over the grid, so nothing is downloaded.
runMulti <- function(out, reps, end = 4) {
  polys <- terra::as.polygons(terra::ext(testGrid()), crs = "EPSG:3400")
  polys$YEAR <- 2000L
  testthat::local_mocked_bindings(
    fetch_nbac_polys = function(...) polys,
    fetch_nfdb_polys = function(...) polys,
    .package = "fireregimetools"
  )
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = end),
    params = list(burnSummaries = list(
      mode = "multi", reps = reps, simTimes = c(0, end), summaryPeriod = c(2L, 4L),
      summaryInterval = 1L, .studyAreaName = "test", .plots = "png"
    )),
    modules = "burnSummaries",
    objects = list(reportingPolygons = polys),
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )
}

test_that("multi mode with every replicate fire-free keeps the fire-size columns", {
  mods <- writeLandMineStandIn(withr::local_tempdir())

  ## the columns and types multi mode reads back from a replicate with fires
  withFire <- file.path(testPaths$outputPath, "multi-withFire")
  runSingleRep(file.path(withFire, "rep01"), rep = 1L, fires = TRUE, mods = mods)
  ref <- data.table::as.data.table(dplyr::collect(
    fireregimetools::open_burn_dataset(file.path(withFire, "rep01", "burnSummaries_fireSizes"))
  ))
  expect_gt(nrow(ref), 0L)

  allFree <- file.path(testPaths$outputPath, "multi-fireFree")
  for (i in 1:2) {
    runSingleRep(file.path(allFree, sprintf("rep%02d", i)), rep = i, fires = FALSE, mods = mods)
  }
  expect_null(fireregimetools::open_burn_dataset(
    file.path(allFree, sprintf("rep%02d", 1:2), "burnSummaries_fireSizes")
  ))

  expect_no_warning(sim <- runMulti(allFree, reps = 1:2), message = "Input has no columns")

  expect_identical(nrow(sim$fireSizes), 0L)
  expect_identical(lapply(sim$fireSizes, class), lapply(ref, class))
  expect_identical(
    readLines(file.path(allFree, "burnSummaries_fireSizes_allReps.csv")),
    paste(names(ref), collapse = ",")
  )

  ## the burn maps are still drawn; the fire-size histogram, with no fires, is not
  figs <- basename(list.files(allFree, "\\.png$", recursive = TRUE))
  expect_true("cumulative_burn_maps.png" %in% figs)
  expect_false("simulated_number_size_fires.png" %in% figs)
})
