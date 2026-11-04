## summaryPeriod defaults to the last 30% of a 1000-year run (years 700-1000). A shorter run has
## no such years, so the default adapts to the run; a supplied period outside the run is an error.

test_that("a 1000-year run keeps the default 700-1000", {
  expect_identical(resolveSummaryPeriod(c(NA_integer_, NA_integer_), start = 0, end = 1000),
                   c(700L, 1000L))
  expect_identical(resolveSummaryPeriod(c(NA_integer_, NA_integer_), start = 0, end = 1500),
                   c(700L, 1000L))
})

test_that("a short run defaults to its last 30%, with a message", {
  expect_message(p <- resolveSummaryPeriod(c(NA_integer_, NA_integer_), start = 0, end = 10),
                 "last 30%")
  expect_identical(p, c(7L, 10L))
  expect_message(p <- resolveSummaryPeriod(c(NA_integer_, NA_integer_), start = 2020, end = 2120),
                 "2090")
  expect_identical(p, c(2090L, 2120L))
})

test_that("a supplied period inside the run is kept, without a message", {
  expect_no_message(p <- resolveSummaryPeriod(c(2L, 8L), start = 0, end = 10))
  expect_identical(p, c(2L, 8L))
})

test_that("a supplied period outside the run is an error naming both ranges", {
  expect_error(resolveSummaryPeriod(c(700L, 1000L), start = 0, end = 10),
               "summaryPeriod.*700.*1000.*0.*10")
})

test_that("init stops on a supplied period outside the run", {
  sim <- suppressMessages(SpaDES.core::simInit(
    times = list(start = 0, end = 10), modules = moduleName, paths = testPaths,
    objects = list(reportingPolygons = 1),
    params = list(burnSummaries = list(mode = "multi", summaryPeriod = c(700L, 1000L)))))
  expect_error(suppressMessages(SpaDES.core::spades(sim, debug = FALSE)),
               "summaryPeriod.*700.*1000.*0.*10")
})

## summaryPeriod is absolute simulation times everywhere: with start = 2020 and summaryPeriod
## c(2020, 3020) the summary times are 2020 to 3020, not 4040 to 5040.
test_that("summary times are the period itself, by summaryInterval", {
  expect_equal(summaryTimes(c(2020L, 3020L), 100L), seq(2020, 3020, by = 100))
})

test_that("single-mode init schedules the period's first save at summaryPeriod[1], not start + it", {
  r <- terra::rast(nrows = 2, ncols = 2, vals = 1L, crs = "EPSG:3978")
  sim <- suppressMessages(SpaDES.core::simInit(
    times = list(start = 2020, end = 3020), modules = moduleName, paths = testPaths,
    objects = list(burnMap = r, flammableMap = r, rstTimeSinceFire = r, rstCurrentBurn = r),
    params = list(burnSummaries = list(mode = "single", summaryPeriod = c(2720L, 3020L)))))
  sim <- suppressMessages(SpaDES.core::spades(sim, debug = FALSE, events = list(burnSummaries = "init")))
  ev <- SpaDES.core::events(sim)
  saves <- ev[ev$eventType == "save_single", ]
  expect_true(2720 %in% saves$eventTime)
  expect_false(any(saves$eventTime > 3020))
})
