## update_tsf runs after the fire model's burn in each year, so a cell that burns in a year reads 0
## that year, as the fire model's own count does (scfmSpread's timeSinceFire). It ran at the default
## priority (5), after LandMine's Burn (2.5) but before fireSense_burn's burn (5.13) and scfmSpread's
## (7.5), so with those it reset a burned cell a year late.

## A stand-in fire model that burns cell 1 in year 2 only, at `priority`.
writeOnceBurner <- function(dir, priority) {
  d <- file.path(dir, "onceBurner")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "defineModule(sim, list(",
      '  name = "onceBurner", description = "test stand-in for a fire model",',
      '  keywords = "test", authors = person("Test", "Author"), childModules = character(0),',
      '  version = list(onceBurner = "0.0.1"), timeframe = as.POSIXlt(c(NA, NA)),',
      '  timeunit = "year", citation = list(), documentation = list(),',
      '  reqdPkgs = list("data.table", "terra"),',
      "  parameters = rbind(),",
      "  inputObjects = bindrows(",
      '    expectsInput("flammableMap", "SpatRaster", "flammable cells")',
      "  ),",
      "  outputObjects = bindrows(",
      '    createsOutput("burnMap", "SpatRaster", "cumulative burns"),',
      '    createsOutput("fireSizes", "list", "fire sizes by year"),',
      '    createsOutput("rstCurrentBurn", "SpatRaster", "cells burned this year")',
      "  )",
      "))",
      "doEvent.onceBurner <- function(sim, eventTime, eventType) {",
      "  switch(eventType,",
      "    init = {",
      "      sim$rstCurrentBurn <- terra::rast(sim$flammableMap, vals = 0L)",
      "      sim$burnMap <- terra::rast(sim$flammableMap, vals = 0L)",
      "      sim$fireSizes <- list()",
      sprintf('      sim <- scheduleEvent(sim, start(sim), "onceBurner", "burn", %s)', priority),
      "    },",
      "    burn = {",
      "      sim$rstCurrentBurn[] <- 0L",
      "      if (time(sim) == 2) {",
      "        sim$rstCurrentBurn[1] <- 1L",
      "        sim$burnMap <- sim$burnMap + sim$rstCurrentBurn",
      '        sim$fireSizes[["2"]] <- data.table::data.table(size = 1L, maxSize = 1L)',
      "      }",
      sprintf('      sim <- scheduleEvent(sim, time(sim) + 1, "onceBurner", "burn", %s)', priority),
      "    }",
      "  )",
      "  invisible(sim)",
      "}"
    ),
    file.path(d, "onceBurner.R")
  )
  dir
}

## Cells 1-3 flammable, starting time since fire 10, 20, 30; cell 4 not flammable (NA). Saves the
## map every year from 1 to 4.
runOnceBurner <- function(out, priority) {
  mods <- writeOnceBurner(withr::local_tempdir(.local_envir = parent.frame()), priority)
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 500, ymin = 0, ymax = 500, crs = "EPSG:3400")
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 4),
    params = list(burnSummaries = list(
      mode = "single", reps = 1L, summaryPeriod = c(1L, 4L),
      summaryInterval = 1L, .studyAreaName = "test"
    )),
    modules = c("onceBurner", "burnSummaries"),
    loadOrder = c("onceBurner", "burnSummaries"),
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

savedTSFat <- function(out, year) {
  terra::values(terra::rast(file.path(out, paste0("rstTimeSinceFire_year", year, ".tif"))), mat = FALSE)
}

## the fire models' burn priorities: LandMine's Burn, fireSense_burn's burn, scfmSpread's burn
for (priority in c(LandMine = 2.5, fireSense_burn = 5.13, scfmSpread = 7.5)) {
  test_that(paste("a cell burned in a year reads 0 that year, burning at priority", priority), {
    out <- file.path(testPaths$outputPath, paste0("onceBurner-", priority))
    runOnceBurner(out, priority)

    ## each year adds 1, from year 0; cell 1 burns in year 2
    expect_equal(savedTSFat(out, 1), c(12, 22, 32, NA))
    expect_equal(savedTSFat(out, 2), c(0, 23, 33, NA))
    expect_equal(savedTSFat(out, 3), c(1, 24, 34, NA))
    expect_equal(savedTSFat(out, 4), c(2, 25, 35, NA))
  })
}
