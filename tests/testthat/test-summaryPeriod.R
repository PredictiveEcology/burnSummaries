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
