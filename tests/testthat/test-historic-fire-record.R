## The historic fire record behind the historic burn map (multi mode): NBAC perimeters from NBAC's
## first year (1972) on, NFDB polygons for the years before it. Through 1.1.0.9001 the NFDB filled
## any year with no NBAC fire inside the study area, so a year from 1972 on with no NBAC fire there
## took NFDB polygons, and NFDB polygons with no year (YEAR = -9999) were kept and counted.

## Polygons with the given fire years, shaped as fireregimetools' loaders return them.
yearPolys <- function(years) {
  v <- terra::vect(
    rep("POLYGON ((0 0, 100 0, 100 100, 0 100, 0 0))", length(years)),
    crs = "EPSG:3400"
  )
  v$YEAR <- as.integer(years)
  v$SIZE_HA <- rep(1, length(years))
  v
}

## One polygon inside each of the cells `cells` of raster `r`, with fire years `years`.
cellPolys <- function(r, cells, years) {
  xy <- terra::xyFromCell(r, cells)
  h <- 0.4 * terra::res(r)
  wkt <- sprintf(
    "POLYGON ((%1$f %3$f, %2$f %3$f, %2$f %4$f, %1$f %4$f, %1$f %3$f))",
    xy[, 1] - h[1], xy[, 1] + h[1], xy[, 2] - h[2], xy[, 2] + h[2]
  )
  v <- terra::vect(wkt, crs = terra::crs(r))
  v$YEAR <- as.integer(years)
  v$SIZE_HA <- rep(1, length(years))
  v
}

test_that("NFDB fills only the years before NBAC begins, not later years NBAC lacks here", {
  nbac <- yearPolys(c(1980, 1995))
  nfdb <- yearPolys(c(1950, 1971, 1972, 1985, 1980))
  fires <- suppressMessages(historicFirePolys(nbac, nfdb, nbacStartYear = 1972L))
  expect_identical(names(fires), "YEAR")
  ## NBAC has no fire here in 1972 or 1985, but it covers those years, so the NFDB's are not used
  expect_identical(sort(fires$YEAR), c(1950L, 1971L, 1980L, 1995L))
})

test_that("NFDB polygons with no year are dropped", {
  fires <- suppressMessages(historicFirePolys(
    yearPolys(1990),
    yearPolys(c(-9999, NA, 0, 1960)),
    nbacStartYear = 1972L
  ))
  expect_identical(sort(fires$YEAR), c(1960L, 1990L))
})

test_that("nbacStartYear is the boundary for both records, so each year has one source", {
  fires <- suppressMessages(historicFirePolys(
    yearPolys(c(1975, 1985)),
    yearPolys(c(1975, 1979, 1980)),
    nbacStartYear = 1980L
  ))
  ## 1975 from the NFDB only; 1980 from NBAC, which has no fire then
  expect_identical(sort(fires$YEAR), c(1975L, 1979L, 1985L))
})

test_that("either record may be empty", {
  nbacOnly <- suppressMessages(historicFirePolys(yearPolys(1990), yearPolys(c(1985, -9999)), 1972L))
  expect_identical(nbacOnly$YEAR, 1990L)
  nfdbOnly <- suppressMessages(historicFirePolys(yearPolys(1990)[0, ], yearPolys(1960), 1972L))
  expect_identical(nfdbOnly$YEAR, 1960L)
  expect_identical(names(nfdbOnly), "YEAR")
})

test_that("the log names the backfilled years and the NFDB polygons left out", {
  msgs <- testthat::capture_messages(
    historicFirePolys(yearPolys(1980), yearPolys(c(1971, 1960, 1990, -9999, NA)), 1972L)
  )
  expect_match(msgs, "backfilling 2 NFDB-only year\\(s\\): 1960, 1971[[:space:]]*$", all = FALSE)
  expect_match(msgs, "not used: 2 with no year, 1 from 1972 on", all = FALSE)
})

test_that("the default boundary is NBAC's first year, 1972", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(md$parameters$default[[which(md$parameters$paramName == "nbacStartYear")]], 1972L)
})

test_that("the historic burn map takes NFDB only before 1972 and drops undated polygons", {
  r <- terra::rast(nrows = 1, ncols = 7, xmin = 0, xmax = 700, ymin = 0, ymax = 100, crs = "EPSG:3400")
  out <- file.path(testPaths$outputPath, "multi-historic")
  dir.create(file.path(out, "rep01"), recursive = TRUE, showWarnings = FALSE)
  ## what multi mode reads from each replicate; no simulated fire
  terra::writeRaster(terra::rast(r, vals = 1L), file.path(out, "rep01", "flammableMap_year4.tif"))
  terra::writeRaster(terra::rast(r, vals = 0L), file.path(out, "rep01", "burnMap_year4.tif"))

  ## cell 1: NBAC 1980. NFDB: cell 2 1971, before NBAC begins; cell 3 1972, NBAC's first year, with
  ## no NBAC fire here; cell 4 1990, a later year with no NBAC fire here; cells 5 and 6 no year
  ## (-9999, NA); cell 7 1980, a year NBAC has here.
  nbac <- cellPolys(r, 1, 1980)
  nfdb <- cellPolys(r, 2:7, c(1971, 1972, 1990, -9999, NA, 1980))
  testthat::local_mocked_bindings(
    fetch_nbac_polys = function(...) nbac,
    fetch_nfdb_polys = function(...) nfdb,
    .package = "fireregimetools"
  )
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 4),
    params = list(burnSummaries = list(
      mode = "multi", reps = 1L, simTimes = c(0, 4), summaryPeriod = c(2L, 4L),
      summaryInterval = 1L, .studyAreaName = "test", .plots = NA_character_
    )),
    modules = "burnSummaries",
    objects = list(reportingPolygons = nbac),
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )

  hist <- terra::values(
    terra::rast(file.path(out, "meanAnnualCumulBurnMapHistoric.tif")),
    mat = FALSE
  )
  ## one fire each in cells 1 (NBAC 1980) and 2 (NFDB 1971), over a record of two fire years
  expect_identical(which(hist > 0), 1:2)
  expect_equal(hist[1:2], c(0.5, 0.5))
})
