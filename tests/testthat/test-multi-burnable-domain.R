## Multi mode's two burn maps cover the same cells, and the historic map is a rate per year of record.
## Through 1.1.0.9002 the historic map covered every cell of the grid's bounding box with a recorded
## fire and the simulated map every flammable cell, while LandMine burns only inside its study area
## and ignites only where there is a fire return interval. The historic map was divided by the number
## of distinct years with a fire anywhere in the bounding box, which leaves fire-free years out.

## A 1 x 8 grid of 100 m cells. Cell 8 is not flammable. The study area holds the centres of cells
## 1-6. Cell 3 has no fire return interval and cell 4 an interval of 0. That leaves cells 1, 2, 5 and
## 6 burnable.
domainGrid <- function() {
  terra::rast(nrows = 1, ncols = 8, xmin = 0, xmax = 800, ymin = 0, ymax = 100, crs = "EPSG:3400")
}
domainStudyArea <- function() {
  terra::vect("POLYGON ((0 0, 600 0, 600 100, 0 100, 0 0))", crs = "EPSG:3400")
}
domainFRI <- function() {
  terra::rast(domainGrid(), vals = c(100L, 100L, NA, 0L, 100L, 100L, 100L, 100L))
}
burnable <- c(1L, 2L, 5L, 6L)

## The fire record. In the burnable cells: NFDB 1961 (cell 5), NBAC 1980 and 1990 (cell 1) and 2000
## (cell 2), so the record there runs 1961-2000, 40 years, only 4 of them with a fire. Elsewhere,
## fires that must not stretch it: NBAC 2005 (cell 3, no interval) and 2010 (cell 7, outside the
## study area); NFDB 1940 (cell 7) and 1950 (cell 8, not flammable). In the grid's bounding box the
## record runs 1940-2010, with 8 distinct fire years.
domainFires <- function(r) {
  list(
    nbac = cellPolys(r, c(1, 1, 2, 3, 7), c(1980, 1990, 2000, 2005, 2010)),
    nfdb = cellPolys(r, c(5, 7, 8), c(1961, 1940, 1950))
  )
}

## Multi mode over one 4-year replicate, with the fire records mocked and `objects` added to the
## simulation and `params` to burnSummaries' parameters. Returns the two maps it writes (NULL for a map
## not written), the observed row of the fire regime table and the messages it logs.
runMultiDomain <- function(out, objects = list(), params = list()) {
  r <- domainGrid()
  dir.create(file.path(out, "rep01"), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(
    terra::rast(r, vals = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L)),
    file.path(out, "rep01", "flammableMap_year4.tif")
  )
  terra::writeRaster(
    terra::rast(r, vals = c(4L, 2L, 4L, 4L, 0L, 1L, 4L, 4L)),
    file.path(out, "rep01", "burnMap_year4.tif")
  )
  fires <- domainFires(r)
  testthat::local_mocked_bindings(
    fetch_nbac_polys = function(...) fires$nbac,
    fetch_nfdb_polys = function(...) fires$nfdb,
    .package = "fireregimetools"
  )
  msgs <- testthat::capture_messages(SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 4),
    params = list(burnSummaries = utils::modifyList(list(
      mode = "multi", reps = 1L, simTimes = c(0, 4), summaryPeriod = c(2L, 4L),
      summaryInterval = 1L, .studyAreaName = "test", .plots = NA_character_
    ), params)),
    modules = "burnSummaries",
    objects = objects,
    paths = utils::modifyList(testPaths, list(outputPath = out))
  ))
  readMap <- function(f) {
    f <- file.path(out, f)
    if (file.exists(f)) terra::values(terra::rast(f), mat = FALSE)
  }
  list(
    historic = readMap("meanAnnualCumulBurnMapHistoric.tif"),
    simulated = readMap("meanAnnualCumulBurnMap.tif"),
    observed = subset(data.table::fread(file.path(out, "burnSummaries_fireRegime.csv")), source == "observed"),
    messages = msgs
  )
}

test_that("both burn maps cover only the burnable cells: in the study area, with an interval", {
  maps <- runMultiDomain(
    file.path(testPaths$outputPath, "multi-domain"),
    list(studyArea = domainStudyArea(), fireReturnInterval = domainFRI())
  )
  expect_identical(which(!is.na(maps$historic)), burnable)
  expect_identical(which(!is.na(maps$simulated)), burnable)
  ## burns over 4 years; cell 6 had no recorded fire and cell 5 no simulated one: rates of 0, not NA
  expect_equal(maps$simulated[burnable], c(4, 2, 0, 1) / 4, tolerance = 1e-6)
  expect_identical(maps$historic[6], 0)
})

test_that("the historic rate divides by the record's span in the burnable cells", {
  maps <- runMultiDomain(
    file.path(testPaths$outputPath, "multi-span"),
    list(studyArea = domainStudyArea(), fireReturnInterval = domainFRI())
  )
  ## 1961-2000: 40 years, fire-free ones included. Not the 4 years with a fire in these cells, nor the
  ## 8 in the bounding box, nor its 1940-2010 span.
  expect_equal(maps$historic[burnable], c(2, 1, 1, 0) / 40, tolerance = 1e-6)
  expect_match(maps$messages, "historic record 1961-2000 \\(40 years\\)", all = FALSE)
})

test_that("without studyArea or fireReturnInterval the maps cover every flammable cell, and say so", {
  maps <- runMultiDomain(file.path(testPaths$outputPath, "multi-noDomain"))
  expect_identical(which(!is.na(maps$historic)), 1:7)
  expect_identical(which(!is.na(maps$simulated)), 1:7)
  ## the record in the flammable cells runs 1940-2010 (cell 7); cell 8's 1950 fire is left out
  expect_equal(maps$historic[1:7], c(2, 1, 1, 0, 1, 0, 2) / 71, tolerance = 1e-6)
  expect_equal(maps$simulated[1:7], c(4, 2, 4, 4, 0, 1, 4) / 4, tolerance = 1e-6)
  expect_match(maps$messages, "historic record 1940-2010 \\(71 years\\)", all = FALSE)
  expect_match(maps$messages, "no `studyArea`", all = FALSE)
  expect_match(maps$messages, "no `fireReturnInterval`", all = FALSE)
})

test_that("a study area in another CRS, or as sf, selects the same cells; one off the grid stops", {
  flm <- terra::rast(domainGrid(), vals = c(rep(1L, 7), 0L))
  sa <- domainStudyArea()
  ref <- suppressMessages(burnableCells(flm, sa))
  expect_identical(which(ref), 1:6)
  lonlat <- sf::st_as_sf(terra::project(sa, "EPSG:4326"))
  expect_identical(suppressMessages(burnableCells(flm, lonlat)), ref)
  expect_error(
    suppressMessages(burnableCells(flm, terra::shift(sa, dx = 1e5))),
    "No burnable cell"
  )
})

test_that("historicYears sets the record on the burnable cells: its fires and its span", {
  maps <- runMultiDomain(
    file.path(testPaths$outputPath, "multi-window"),
    list(studyArea = domainStudyArea(), fireReturnInterval = domainFRI()),
    list(historicYears = c(1971L, 2010L))
  )
  ## 1971-2010: 40 years. NFDB 1961 (cell 5) is outside the window; the 2005 and 2010 fires are in it
  ## but not on burnable cells, so they count in neither the map nor the span.
  expect_equal(maps$historic[burnable], c(2, 1, 0, 0) / 40, tolerance = 1e-6)
  expect_match(maps$messages, "historic record 1971-2010 \\(40 years\\)", all = FALSE)
  expect_identical(c(maps$observed$firstYear, maps$observed$lastYear), c(1971L, 2010L))
  r <- domainGrid()
  expect_equal(maps$observed$areaBurnedHaPerYear, 3 * prod(terra::res(r)) / 1e4 / 40, tolerance = 1e-6)
})

test_that("supplied firePolys replace the download, on the burnable cells, without undated fires", {
  r <- domainGrid()
  supplied <- list(cellPolys(r, c(1, 5, 6, 7), c(1990, 1995, -9999, 2000)))
  maps <- runMultiDomain(
    file.path(testPaths$outputPath, "multi-supplied"),
    list(studyArea = domainStudyArea(), fireReturnInterval = domainFRI(), firePolys = supplied)
  )
  ## the record in the burnable cells runs 1990-1995 (6 years); cell 6's undated fire is dropped and
  ## cell 7 is outside the study area
  expect_equal(maps$historic[burnable], c(1, 0, 1, 0) / 6, tolerance = 1e-6)
  expect_identical(which(!is.na(maps$historic)), burnable)
  expect_match(maps$messages, "historic record 1990-1995 \\(6 years\\)", all = FALSE)
  ## nothing was fetched
  expect_false(any(grepl("NBAC perimeters", maps$messages)))
})

test_that("the observed fire count is of fires on the burnable cells", {
  maps <- runMultiDomain(
    file.path(testPaths$outputPath, "multi-nFires"),
    list(studyArea = domainStudyArea(), fireReturnInterval = domainFRI())
  )
  ## 8 recorded fires; NBAC 1980, 1990, 2000 and NFDB 1961 are on burnable cells
  expect_identical(maps$observed$nFires, 4L)
})
