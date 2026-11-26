## Awkward observed fire polygons: a fire mapped as several parts, one straddling the edge of the grid, and
## supplied polygons that are invalid, empty or carry a YEAR that is text.

## A 4 x 4 grid of 1 ha cells (100 m; EPSG:3400 units are metres), every cell burnable.
recordGrid <- function() {
  terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 400, ymin = 0, ymax = 400, crs = "EPSG:3400", vals = 1L)
}

test_that("a fire mapped as one multipart polygon is one fire, its size all its parts", {
  r <- recordGrid()
  ## squares around the centres of cell 1 (top left) and cell 16 (bottom right)
  p <- terra::vect(
    "MULTIPOLYGON (((10 310, 90 310, 90 390, 10 390, 10 310)), ((310 10, 390 10, 390 90, 310 90, 310 10)))",
    crs = "EPSG:3400"
  )
  p$YEAR <- 2000L
  yrs <- 2000:2009
  h <- historicBurnMap(p, r, rep(TRUE, 16), yrs)
  expect_equal(which(terra::values(h, mat = FALSE) > 0), c(1L, 16L))
  o <- observedFireRegime(p, h, r, yrs)
  expect_identical(o$nFires, 1L)
  expect_equal(o$largestFireHa, 2)
  expect_equal(o$areaBurnedHaPerYear, 2 / 10)
})

test_that("a fire straddling the edge of the grid counts only the land inside it", {
  r <- recordGrid()
  ## from the centre of cell 4 (top right) to 100 m beyond the grid
  p <- terra::vect("POLYGON ((310 310, 500 310, 500 390, 310 390, 310 310))", crs = "EPSG:3400")
  p$YEAR <- 1995L
  yrs <- suppressMessages(historicRecordYears(p))
  expect_identical(yrs, 1995L)
  h <- historicBurnMap(p, r, rep(TRUE, 16), yrs)
  expect_equal(which(terra::values(h, mat = FALSE) > 0), 4L)
  o <- observedFireRegime(p, h, r, yrs)
  expect_identical(o$nFires, 1L)
  expect_equal(o$largestFireHa, 1)
})

test_that("supplied fire polygons are repaired, and a YEAR given as text is read as a number", {
  mods <- writeLandMineStandIn(withr::local_tempdir())
  out <- file.path(testPaths$outputPath, "multi-supplied-awkward")
  runSingleRep(file.path(out, "rep01"), rep = 1L, fires = FALSE, mods = mods)
  ## on testGrid(): 2 x 2 cells of 250 m, cells 1-3 flammable
  polys <- terra::vect(c(
    "POLYGON ((10 260, 240 260, 240 490, 10 490, 10 260))", ## cell 1, YEAR "2000"
    "POLYGON ((0 500, 250 250, 500 500, 250 250, 0 500))", ## invalid and of no area, "2001"
    "POLYGON EMPTY", ## "2002"
    "POLYGON ((260 260, 490 260, 490 490, 260 490, 260 260))" ## cell 2, YEAR "unknown"
  ), crs = "EPSG:3400")
  polys$YEAR <- c("2000", "2001", "2002", "unknown")
  msgs <- testthat::capture_messages(
    runMulti(out, reps = 1L, plots = NA_character_, objects = list(firePolys = list(polys)))
  )
  fr <- data.table::fread(file.path(out, "burnSummaries_fireRegime.csv"))
  obs <- fr[fr$source == "observed", ]
  ## only the 2000 fire: the polygon of no area and the empty one are dropped by the repair, and the one with
  ## no year is left out. Unrepaired, the polygon of no area had no cell in the map but counted as a fire.
  expect_identical(obs$nFires, 1L)
  expect_identical(c(obs$firstYear, obs$lastYear), c(2000L, 2000L))
  expect_match(msgs, "1 fire polygon\\(s\\) with no year", all = FALSE)
  h <- terra::values(terra::rast(file.path(out, "meanAnnualCumulBurnMapHistoric.tif")), mat = FALSE)
  expect_equal(h[1], 1) ## cell 1 burned in the one year of record
  expect_equal(h[2], 0)
})
