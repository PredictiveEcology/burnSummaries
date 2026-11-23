## 10 x 10 grid of 240 m cells (EPSG:3400 units are metres); minCellHa 100 -> factor 5
fineGrid <- function(n = 10, vals = 1) {
  terra::rast(nrows = n, ncols = n, vals = vals, crs = "EPSG:3400",
              xmin = 0, xmax = n * 240, ymin = 0, ymax = n * 240)
}
polyOf <- function(r, cells, years, cause = "Natural", prescribed = "false") {
  xy <- terra::xyFromCell(r, cells)
  h <- 0.4 * terra::res(r)
  wkt <- sprintf("POLYGON ((%1$f %3$f, %2$f %3$f, %2$f %4$f, %1$f %4$f, %1$f %3$f))",
                 xy[, 1] - h[1], xy[, 1] + h[1], xy[, 2] - h[2], xy[, 2] + h[2])
  v <- terra::vect(wkt, crs = terra::crs(r))
  v$YEAR <- as.integer(years)
  v$FIRECAUS <- cause
  v$PRESCRIBED <- prescribed
  v
}

test_that("the aggregation factor is the smallest integer whose cell is at least minCellHa", {
  expect_identical(cellFactor(c(240, 240), 100), 5L) ## 4 -> 23 ha, 5 -> 144 ha
  expect_identical(cellFactor(c(100, 100), 100), 10L) ## exactly 100 ha
  expect_identical(cellFactor(c(2000, 2000), 100), 1L) ## already coarser
})

test_that("prescribed burns and causes outside fireCauses are dropped, with a message", {
  r <- fineGrid()
  p <- polyOf(r, 1:4, 2000L, cause = c("Natural", "Undetermined", "Human", "Natural"),
              prescribed = c("false", "false", "false", "true"))
  kept <- fireSenseUtils::modelledFires(p, fireSenseUtils::fireSenseFireCauses, "nbac")
  expect_equal(nrow(kept), 2) ## natural and undetermined; the human fire and the prescribed one go
  ## the human cause is kept when the fits model it
  kept2 <- fireSenseUtils::modelledFires(p, c("natural", "unknown", "human"), "nbac")
  expect_equal(nrow(kept2), 3) ## still no prescribed
})

test_that("historicFirePolys keeps only modelled-cause, non-prescribed fires of both records", {
  r <- fineGrid()
  nbac <- polyOf(r, 1:3, 1990L, cause = c("Natural", "Human", "Undetermined"),
                 prescribed = c("false", "false", "true"))
  nfdb <- polyOf(r, 4:6, 1960L, prescribed = c(NA, NA, "yes")) ## NFDB polygons: PRESCRIBED "yes"
  nfdb$CAUSE <- c("N", "H", "N")
  fires <- suppressMessages(historicFirePolys(nbac, nfdb, 1972L, fireSenseUtils::fireSenseFireCauses))
  expect_identical(sort(fires$YEAR), c(1960L, 1990L)) ## NBAC cell 1, NFDB cell 4; the prescribed NFDB natural fire (cell 6) goes
})

test_that("the historic map is the share of burnable area burned per year, two fires on one patch count once", {
  r <- fineGrid()
  ## coarse cell 1 = fine cells in rows 1-5, cols 1-5. Year 2000: two overlapping polygons on fine cell 1
  ## (a second one on the same cell), year 2001: fine cells 1 and 2.
  p <- rbind(polyOf(r, c(1, 1, 1, 2), c(2000L, 2000L, 2001L, 2001L)))
  burnable <- rep(TRUE, terra::ncell(r))
  yrs <- 2000:2009
  mm <- historicBurnMap(p, r, burnable, yrs, fact = 5L)
  expect_equal(dim(mm)[1:2], c(2, 2))
  ## year 2000: 1 of 25 cells; 2001: 2 of 25; over 10 years
  expect_equal(terra::values(mm, mat = FALSE)[1], (1 / 25 + 2 / 25) / 10)
  expect_equal(terra::values(mm, mat = FALSE)[2:4], c(0, 0, 0))
})

test_that("coarse cells with no burnable fine cell are NA, and burnable area is the denominator", {
  r <- fineGrid()
  burnable <- rep(FALSE, terra::ncell(r))
  rc <- terra::rowColFromCell(r, seq_len(terra::ncell(r)))
  burnable[rc[, 1] <= 5 & rc[, 2] <= 5 & (rc[, 2] <= 2)] <- TRUE ## 10 burnable cells in coarse cell 1
  p <- polyOf(r, 1, 2000L)
  mm <- historicBurnMap(p, r, burnable, 2000:2009, fact = 5L)
  v <- terra::values(mm, mat = FALSE)
  expect_equal(v[1], (1 / 10) / 10)
  expect_true(all(is.na(v[2:4])))
})

test_that("the simulated map is the mean of the burnable fine-cell rates in each coarse cell", {
  r <- fineGrid()
  burnable <- rep(TRUE, terra::ncell(r))
  burnable[2] <- FALSE
  rate <- terra::rast(r, vals = 0)
  rate[1] <- 0.5; rate[3] <- 0.25
  sim <- coarseShare(rate, burnable, 5L)
  expect_equal(terra::values(sim, mat = FALSE)[1], (0.5 + 0.25) / 24)
})

test_that("multi mode writes both maps at the coarse cell and drops non-modelled fires", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2400, ymin = 0, ymax = 2400, crs = "EPSG:3400")
  out <- file.path(testPaths$outputPath, "multi-coarse")
  dir.create(file.path(out, "rep01"), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(terra::rast(r, vals = 1L), file.path(out, "rep01", "flammableMap_year4.tif"), overwrite = TRUE)
  terra::writeRaster(terra::rast(r, vals = 0L), file.path(out, "rep01", "burnMap_year4.tif"), overwrite = TRUE)
  nbac <- polyOf(r, c(1, 1, 1, 2), c(2000L, 2000L, 2001L, 2002L),
                 cause = c("Natural", "Undetermined", "Natural", "Human"),
                 prescribed = c("false", "false", "false", "false"))
  nbac$PRESCRIBED[3] <- "true"
  nbac$CAUSE <- "N"
  nfdb <- nbac[0, ]
  testthat::local_mocked_bindings(
    fetch_nbac_polys = function(...) nbac,
    fetch_nfdb_polys = function(...) nfdb,
    .package = "fireregimetools"
  )
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 4),
    params = list(burnSummaries = list(
      mode = "multi", reps = 1L, simTimes = c(0, 4), summaryPeriod = c(2L, 4L), summaryInterval = 1L,
      historicYears = c(2000L, 2009L), .studyAreaName = "test", .plots = NA_character_
    )),
    modules = "burnSummaries",
    paths = utils::modifyList(testPaths, list(outputPath = out))
  )
  h <- terra::rast(file.path(out, "meanAnnualCumulBurnMapHistoric.tif"))
  s <- terra::rast(file.path(out, "meanAnnualCumulBurnMap.tif"))
  expect_equal(dim(h)[1:2], c(2, 2))
  expect_equal(dim(s)[1:2], c(2, 2))
  ## only year 2000 (two polygons, one cell) survives: 1 of 25 cells, over 10 years
  expect_equal(terra::values(h, mat = FALSE)[1], (1 / 25) / 10, tolerance = 1e-6)
  expect_equal(terra::values(s, mat = FALSE), rep(0, 4))
})

test_that("a record with no cause column is an error, as in the fit", {
  p <- polyOf(fineGrid(), 1:2, 2000L)
  p$FIRECAUS <- NULL
  expect_error(historicFirePolys(p, p[integer(0), ], 1972L, "natural"), "FIRECAUS")
})
