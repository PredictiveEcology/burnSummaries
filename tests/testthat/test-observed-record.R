## Awkward observed fire polygons: a fire mapped as several parts or split at a boundary into several
## polygons, one straddling the edge of the grid, and supplied polygons that are invalid, empty or carry a
## YEAR that is text.

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

test_that("NBAC polygons sharing an NFIREID in a year are one fire; each NFDB polygon is one fire", {
  r <- recordGrid()
  ## NBAC splits a fire at a provincial, territorial or park boundary into polygons that share its
  ## NFIREID in that year. NFIREID 7 in 2000 is one fire in cells 1 and 2; NFIREID 7 in 2001 is another
  ## fire (the ID is assigned per year). Through 2.0.0 each polygon counted as a fire.
  nbac <- cellPolys(r, c(1, 2, 3), c(2000, 2000, 2001))
  nbac$NFIREID <- c(7, 7, 7)
  ## two NFDB fires in one year, which have no NFIREID: two fires, not grouped by their missing ID
  nfdb <- cellPolys(r, c(5, 6), c(1960, 1960))
  fires <- suppressMessages(historicFirePolys(nbac, nfdb, nbacStartYear = 1972L))
  expect_identical(sort(fires$NFIREID, na.last = TRUE), c(7, 7, 7, NA, NA))
  yrs <- 1960:2009
  h <- historicBurnMap(fires, r, rep(TRUE, 16), yrs)
  expect_equal(which(terra::values(h, mat = FALSE) > 0), c(1L, 2L, 3L, 5L, 6L)) ## the map is per cell, as before
  msgs <- testthat::capture_messages(o <- observedFireRegime(fires, h, r, yrs))
  expect_identical(o$nFires, 4L)
  expect_equal(o$largestFireHa, 2) ## 1 ha cells: both pieces of the 2000 fire
  expect_equal(o$areaBurnedHaPerYear, 5 / 50)
  expect_match(msgs, "5 polygons are 4 fires", all = FALSE)

  ## without the IDs, each polygon is a fire
  o1 <- observedFireRegime(fires[, "YEAR"], h, r, yrs)
  expect_identical(o1$nFires, 5L)
  expect_equal(o1$largestFireHa, 1)
})
