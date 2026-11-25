## Checks that turn a silent loss into an error or a warning: a replicate's fire sizes counted for another
## replicate, a fire whose year is not a number, areas on a grid that is not in metres, and fire sizes that
## do not add up to the burns on the burn map.

test_that("the replicate comes from `reps` when it names one, else the last repNN in the path", {
  expect_identical(repIDFromRun(3L, "/out/rep07"), 3L)
  expect_identical(repIDFromRun(1:10, "/out/sa/mainSim/rep07"), 7L)
  expect_identical(repIDFromRun(1:10, "/out/rep07/sub"), 7L) ## NA through 1.1.0.9010
  expect_identical(repIDFromRun(1:10, "/out/rep2023/rep07"), 7L)
  expect_identical(repIDFromRun(1:10, "C:\\out\\rep12"), 12L)
  expect_identical(repIDFromRun(1:10, "/out/representative"), NA_integer_)
})

test_that("each replicate's fire sizes take the replicate of the directory they were read from", {
  expect_identical(repFromDir(data.frame(rep = c(NA, NA), year = c("1", "2"), simSize = 1:2), 3L, "x")$rep,
                   c(3L, 3L))
  expect_identical(repFromDir(data.frame(year = "1", simSize = 1L), 3L, "x")$rep, 3L)
  expect_identical(repFromDir(data.frame(rep = 3L, year = "1", simSize = 1L), 3L, "x")$rep, 3L)
  expect_error(
    repFromDir(data.frame(rep = c(3L, 1L), year = c("1", "2"), simSize = 1:2), 3L, "rep03/fs.csv"),
    "rep03/fs.csv are for replicate 1, not 3"
  )
})

test_that("a fire-size year that is missing or not a number stops", {
  expect_true(checkFireYears(c("1", "25", "1000")))
  expect_true(checkFireYears(c(1, 2)))
  expect_true(checkFireYears(character(0)))
  expect_error(checkFireYears(c("1", "year2")), "year2")
  expect_error(checkFireYears(c(1, NA)), "missing or not a number")
})

test_that("areas need a projected grid in metres", {
  inM <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 200, ymin = 0, ymax = 200, crs = "EPSG:3400")
  expect_silent(checkMetreGrid(inM))
  expect_error(checkMetreGrid(terra::rast(nrows = 2, ncols = 2, crs = "EPSG:4326")), "lon/lat")
  expect_error(checkMetreGrid(terra::rast(nrows = 2, ncols = 2, crs = "EPSG:2263")), "units of 0.3048")
  expect_error(
    checkMetreGrid(terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 200, ymin = 0, ymax = 200, crs = "")),
    "is missing"
  )
})

test_that("each replicate's fire sizes must add up to the burns on its burn map", {
  fs <- data.frame(rep = c(1L, 1L, 2L), simSize = c(3, 2, 4))
  expect_no_warning(checkBurnTotals(fs, c(`1` = 5, `2` = 4)))
  expect_warning(out <- checkBurnTotals(fs, c(`1` = 5, `2` = 6)), class = "burnSummariesBurnTotals")
  expect_equal(out$fireSizes, c(5, 4))
  expect_equal(out$burnMap, c(5, 6))
  ## burns on a replicate's map, and no fire sizes for it
  expect_warning(checkBurnTotals(fs, c(`1` = 5, `2` = 4, `3` = 2)), "3 \\(0 from fire sizes, 2 on the burn map\\)")
})
