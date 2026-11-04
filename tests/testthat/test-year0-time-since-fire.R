## Year 0's time since fire is the map the run starts from. Through 1.1.0 it was saved at .last(),
## after update_tsf had added fireTimestep at start(sim), so the saved year-0 map was the starting
## map + 1. update_tsf itself still runs at start(sim): a fire model reads the starting map + 1 in
## year 1, as before.

## A stand-in fire model: from start(sim) + 1, at LandMine's priority for Burn (2.5), it records the
## time since fire it reads, then burns cell 1.
writeFireStandIn <- function(dir) {
  d <- file.path(dir, "fireStandIn")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "defineModule(sim, list(",
      '  name = "fireStandIn", description = "test stand-in for a fire model",',
      '  keywords = "test", authors = person("Test", "Author"), childModules = character(0),',
      '  version = list(fireStandIn = "0.0.1"), timeframe = as.POSIXlt(c(NA, NA)),',
      '  timeunit = "year", citation = list(), documentation = list(),',
      '  reqdPkgs = list("data.table", "terra"),',
      "  parameters = rbind(),",
      "  inputObjects = bindrows(",
      '    expectsInput("flammableMap", "SpatRaster", "flammable cells"),',
      '    expectsInput("rstTimeSinceFire", "SpatRaster", "time since fire")',
      "  ),",
      "  outputObjects = bindrows(",
      '    createsOutput("burnMap", "SpatRaster", "cumulative burns"),',
      '    createsOutput("fireSizes", "list", "fire sizes by year"),',
      '    createsOutput("rstCurrentBurn", "SpatRaster", "cells burned this year"),',
      '    createsOutput("tsfRead", "list", "time since fire read before each burn, by year")',
      "  )",
      "))",
      "doEvent.fireStandIn <- function(sim, eventTime, eventType) {",
      "  switch(eventType,",
      "    init = {",
      "      sim$rstCurrentBurn <- terra::rast(sim$flammableMap, vals = 0L)",
      "      sim$burnMap <- terra::rast(sim$flammableMap, vals = 0L)",
      "      sim$fireSizes <- list()",
      "      sim$tsfRead <- list()",
      '      sim <- scheduleEvent(sim, start(sim) + 1, "fireStandIn", "burn", 2.5)',
      "    },",
      "    burn = {",
      "      yr <- as.character(time(sim))",
      "      sim$tsfRead[[yr]] <- terra::values(sim$rstTimeSinceFire, mat = FALSE)",
      "      sim$rstCurrentBurn[] <- 0L",
      "      sim$rstCurrentBurn[1] <- 1L",
      "      sim$burnMap <- sim$burnMap + sim$rstCurrentBurn",
      "      sim$fireSizes[[yr]] <- data.table::data.table(size = 1L, maxSize = 1L)",
      '      sim <- scheduleEvent(sim, time(sim) + 1, "fireStandIn", "burn", 2.5)',
      "    }",
      "  )",
      "  invisible(sim)",
      "}"
    ),
    file.path(d, "fireStandIn.R")
  )
  dir
}

## Cells 1-3 flammable, starting time since fire 10, 20, 30; cell 4 not flammable (NA).
tsf0 <- c(10L, 20L, 30L, NA)

runSingle <- function(out, summaryPeriod, summaryInterval, end = 4) {
  mods <- writeFireStandIn(withr::local_tempdir(.local_envir = parent.frame()))
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = "EPSG:3400")
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = end),
    params = list(burnSummaries = list(
      mode = "single", reps = 1L, summaryPeriod = summaryPeriod,
      summaryInterval = summaryInterval, .studyAreaName = "test"
    )),
    modules = c("fireStandIn", "burnSummaries"),
    loadOrder = c("fireStandIn", "burnSummaries"),
    objects = list(
      flammableMap = terra::rast(r, vals = c(1L, 1L, 1L, 0L)),
      rstTimeSinceFire = terra::rast(r, vals = tsf0)
    ),
    paths = utils::modifyList(
      testPaths,
      list(modulePath = c(testPaths$modulePath, mods), outputPath = out)
    )
  )
}

savedTSF <- function(out, year) {
  terra::values(terra::rast(file.path(out, paste0("rstTimeSinceFire_year", year, ".tif"))), mat = FALSE)
}

year0Saves <- function(sim) {
  cmp <- as.data.frame(SpaDES.core::completed(sim))
  sum(cmp$moduleName == "burnSummaries" & cmp$eventType == "save_single" &
    as.numeric(cmp$eventTime) == 0)
}

test_that("year 0 is the starting time since fire, saved once", {
  out <- file.path(testPaths$outputPath, "year0")
  sim <- runSingle(out, summaryPeriod = c(2L, 4L), summaryInterval = 1L)

  expect_equal(savedTSF(out, 0), tsf0)
  expect_identical(year0Saves(sim), 1L)

  ## the fire model's first burn reads the starting map + 1, as before
  expect_equal(sim$tsfRead[["1"]], tsf0 + 1L)

  ## later years unchanged: burned cells reset to 0, the others age by 1 a year
  expect_equal(sim$tsfRead[["2"]], c(0, 22, 32, NA))
  expect_equal(savedTSF(out, 4), c(0, 25, 35, NA))
  expect_false(file.exists(file.path(out, "rstTimeSinceFire_year1.tif")))
  for (year in 2:4) {
    expect_true(file.exists(file.path(out, paste0("rstTimeSinceFire_year", year, ".tif"))))
  }
  expect_equal(terra::values(terra::rast(file.path(out, "burnMap_year4.tif")), mat = FALSE), c(4, 0, 0, 0))
  expect_true(file.exists(file.path(out, "flammableMap_year4.tif")))
})

test_that("a summary period starting at year 0 keeps the starting map and its snapshots", {
  out <- file.path(testPaths$outputPath, "year0-summaryPeriod")
  sim <- runSingle(out, summaryPeriod = c(0L, 4L), summaryInterval = 2L)

  expect_equal(savedTSF(out, 0), tsf0)
  expect_identical(year0Saves(sim), 1L)
  expect_equal(savedTSF(out, 2), c(0, 23, 33, NA))
  expect_equal(savedTSF(out, 4), c(0, 25, 35, NA))
})
