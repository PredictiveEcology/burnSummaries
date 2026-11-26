## `init` in mode "single" must queue each event once per distinct time: `save_single` is wanted
## at the start, at the start of the summary period and at the end, and those times can coincide.

singleSim <- function(summaryPeriod, end = 10, interval = 5L) {
  r <- terra::rast(nrows = 4, ncols = 4, vals = 1L)
  suppressMessages(SpaDES.core::simInit(
    times = list(start = 0, end = end),
    modules = moduleName,
    paths = testPaths,
    params = stats::setNames(
      list(list(mode = "single", summaryPeriod = summaryPeriod, summaryInterval = interval,
                .plots = NA, .useCache = FALSE)),
      moduleName),
    objects = list(burnMap = r, flammableMap = r, rstCurrentBurn = r, rstTimeSinceFire = r)
  ))
}

queuedAfterInit <- function(summaryPeriod, ...) {
  sim <- suppressMessages(SpaDES.core::spades(
    singleSim(summaryPeriod, ...), events = stats::setNames(list("init"), moduleName)))
  ev <- SpaDES.core::events(sim)
  ev[ev$moduleName == moduleName, c("eventTime", "eventType")]
}

test_that("save_single is queued once per distinct time when the summary period starts at the start", {
  ev <- queuedAfterInit(c(0L, 10L))
  saves <- ev$eventTime[ev$eventType == "save_single"]
  expect_identical(sort(saves), c(0, 10))
})

test_that("save_single is queued once per distinct time when the summary period starts at the end", {
  ev <- queuedAfterInit(c(10L, 10L))
  saves <- ev$eventTime[ev$eventType == "save_single"]
  expect_identical(sort(saves), c(0, 10))
})

test_that("save_single keeps its three times when they are distinct", {
  ev <- queuedAfterInit(c(5L, 10L))
  saves <- ev$eventTime[ev$eventType == "save_single"]
  expect_identical(sort(saves), c(0, 5, 10))
})

test_that("save_single runs once per time, with the summary period over the whole simulation", {
  sim <- suppressMessages(SpaDES.core::spades(
    singleSim(c(0L, 10L)), events = stats::setNames(list(c("init", "save_single")), moduleName)))
  done <- SpaDES.core::completed(sim)
  saves <- done$eventTime[done$moduleName == moduleName & done$eventType == "save_single"]
  expect_false(anyDuplicated(saves) > 0L)
  expect_identical(sort(saves), c(0, 5, 10))
})

test_that("simOutputPath is not in the init event's cache key: a run that differs only in it is a cache hit", {
  ## a cache built on one machine is then a hit on another, where the folder differs
  paths <- testPaths
  paths$cachePath <- withr::local_tempdir()
  r <- terra::rast(nrows = 4, ncols = 4, vals = 1L)
  cachedInit <- function(simOutputPath) {
    sim <- suppressMessages(SpaDES.core::simInit(
      times = list(start = 0, end = 10), modules = moduleName, paths = paths,
      params = stats::setNames(list(list(mode = "single", summaryPeriod = c(0L, 10L), .plots = NA,
                                         .useCache = "init", simOutputPath = simOutputPath)), moduleName),
      objects = list(burnMap = r, flammableMap = r, rstCurrentBurn = r, rstTimeSinceFire = r)))
    suppressMessages(SpaDES.core::spades(sim, events = stats::setNames(list("init"), moduleName)))
  }
  nCached <- function() length(unique(reproducible::showCache(paths$cachePath, verbose = -2)$cacheId))
  cachedInit(file.path(tempdir(), "machineA"))
  expect_identical(nCached(), 1L)                                        # the init event, cached
  cachedInit(file.path(tempdir(), "machineB"))
  expect_identical(nCached(), 1L)                                        # a hit: nothing new was cached
})
