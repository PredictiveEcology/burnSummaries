toyRaster <- function(vals) {
  r <- terra::rast(nrows = 2, ncols = 2, vals = vals, crs = "EPSG:3857",
                   xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  r
}

test_that("historic and simulated burn maps share one fill scale", {
  historic <- toyRaster(c(0, 0.1, 0.2, NA))
  simulated <- toyRaster(c(0, 0.5, 1, 0.3))
  p <- cumulBurnMapPlots(historic, simulated, "toy")
  limits <- lapply(1:2, function(i) p[[i]]$scales$get_scales("fill")$limits)
  expect_identical(limits[[1]], limits[[2]])
  expect_equal(limits[[1]], c(0, 1))
})

test_that("historic mean annual burn divides by the record length, not the years with fires", {
  flammable <- toyRaster(rep(1, 4))
  ## fires in 2 of 10 years (2000 and 2009), both covering the first cell
  polys <- terra::vect(
    c("POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))", "POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))"),
    crs = "EPSG:3857"
  )
  polys$YEAR <- c(2000L, 2009L)
  burnable <- rep(TRUE, 4)
  yrs <- suppressMessages(historicRecordYears(polys, flammable, burnable))
  expect_identical(yrs, 2000:2009)
  mm <- historicBurnMap(polys, flammable, burnable, yrs)
  ## a burnable cell with no recorded fire has a rate of 0
  expect_equal(terra::values(mm, mat = FALSE), c(2 / 10, 0, 0, 0))
})

test_that("historicYears sets the record length, so a record with fireless end years is not shortened", {
  flammable <- toyRaster(rep(1, 4))
  polys <- terra::vect("POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))", crs = "EPSG:3857")
  polys$YEAR <- 2000L
  burnable <- rep(TRUE, 4)
  yrs <- suppressMessages(historicRecordYears(polys, flammable, burnable, c(1991L, 2010L)))
  expect_identical(yrs, 1991:2010)
  mm <- historicBurnMap(polys, flammable, burnable, yrs)
  expect_equal(max(terra::values(mm), na.rm = TRUE), 1 / 20)
})

test_that("the fire regime counts each replicate, including one with no fire", {
  fs <- data.frame(rep = c(1, 1, 1, 2), year = c(5, 15, 15, 30), simSize = c(2, 10, 0, 4))
  r <- simulatedFireRegime(fs, reps = 1:3, simTimes = c(0, 40), pixelHa = 5)
  expect_identical(r$rep, 1:3)
  expect_identical(r$nFires, c(2L, 1L, 0L)) ## a fire of 0 cells is not a fire
  expect_equal(r$areaBurnedHaPerYear, c(12 * 5 / 40, 4 * 5 / 40, 0))
  expect_equal(r$largestFireHa, c(50, 20, 0))
})

test_that("the observed regime counts burnable land only", {
  flammable <- toyRaster(c(1, 1, NA, 1)) ## cell 3 cannot burn
  polys <- terra::vect(c("POLYGON ((0 0, 2 0, 2 1, 0 1, 0 0))", "POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))"),
                       crs = "EPSG:3857") ## the bottom row (cells 3, 4) and cell 1
  polys$YEAR <- c(2001L, 2002L)
  yrs <- 2001:2010
  historic <- historicBurnMap(polys, flammable, !is.na(terra::values(flammable, mat = FALSE)), yrs)
  o <- observedFireRegime(polys, historic, flammable, yrs)
  pixHa <- prod(terra::res(flammable)) / 1e4
  expect_equal(o$areaBurnedHaPerYear, 2 * pixHa / 10) ## cells 1 and 4, once each, over 10 years
  expect_equal(o$largestFireHa, pixHa)                ## each fire covers one flammable cell
  expect_identical(c(o$firstYear, o$lastYear), c(2001L, 2010L))
})

test_that("the observed regime has no first or last year without a record", {
  flammable <- toyRaster(rep(1, 4))
  polys <- terra::vect("POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))", crs = "EPSG:3857")[0, ]
  historic <- historicBurnMap(polys, flammable, rep(TRUE, 4), integer(0))
  expect_equal(terra::values(historic, mat = FALSE), rep(0, 4))
  o <- observedFireRegime(polys, historic, flammable, integer(0))
  expect_identical(c(o$firstYear, o$lastYear), c(NA_integer_, NA_integer_))
  expect_identical(o$areaBurnedHaPerYear, 0)
})

test_that("fire per block keeps blocks with no fire as zero", {
  fs <- data.frame(rep = c(1, 1, 2), year = c(5, 25, 6), simSize = c(1, 3, 2))
  b <- fireBlockTable(fs, reps = 1:2, simTimes = c(0, 30), blockYears = 10, pixelHa = 2)
  expect_equal(nrow(b), 2 * 3)
  expect_equal(b$nFires[b$rep == 1], c(1, 0, 1))
  expect_equal(b$areaBurnedHa[b$rep == 1], c(2, 0, 6))
  expect_equal(b$areaBurnedHa[b$rep == 2], c(4, 0, 0))
})

test_that("the simulated regime covers only the years of `period`, default the whole run", {
  fs <- data.frame(rep = c(1, 1, 1, 1), year = c(5, 15, 25, 40), simSize = c(1, 2, 4, 8))
  whole <- simulatedFireRegime(fs, reps = 1, simTimes = c(0, 40), pixelHa = 1)
  expect_equal(whole$nFires, 4L)
  expect_equal(whole$areaBurnedHaPerYear, 15 / 40)
  ## years 10 < year <= 30: the fires of years 15 and 25
  win <- simulatedFireRegime(fs, reps = 1, simTimes = c(0, 40), pixelHa = 1, period = c(10, 30))
  expect_equal(win$nFires, 2L)
  expect_equal(win$areaBurnedHaPerYear, 6 / 20)
  expect_equal(win$largestFireHa, 4)
  expect_identical(c(win$firstYear, win$lastYear), c(10, 30))
})

test_that("the observed regime does not count fires outside the burnable cells", {
  flammable <- toyRaster(c(1, 1, NA, NA)) ## the bottom row cannot burn
  polys <- terra::vect(c("POLYGON ((0 0, 2 0, 2 1, 0 1, 0 0))", "POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))"),
                       crs = "EPSG:3857") ## the bottom row only, and cell 1
  polys$YEAR <- c(2001L, 2002L)
  yrs <- 2001:2010
  historic <- historicBurnMap(polys, flammable, !is.na(terra::values(flammable, mat = FALSE)), yrs)
  o <- observedFireRegime(polys, historic, flammable, yrs)
  expect_identical(o$nFires, 1L)
})
