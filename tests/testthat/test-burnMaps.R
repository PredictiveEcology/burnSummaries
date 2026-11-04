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
  cumul <- terra::rasterize(polys, flammable, field = "YEAR", fun = "count")
  mm <- historicMeanAnnualBurnMap(polys, flammable)
  expect_equal(terra::values(mm, mat = FALSE), terra::values(cumul, mat = FALSE) / 10)
  expect_equal(max(terra::values(mm), na.rm = TRUE), 2 / 10)
})

test_that("historicYears sets the record length, so a record with fireless end years is not shortened", {
  flammable <- toyRaster(rep(1, 4))
  polys <- terra::vect("POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))", crs = "EPSG:3857")
  polys$YEAR <- 2000L
  yrs <- historicRecordYears(polys$YEAR, c(1991L, 2010L))
  expect_identical(yrs, 1991:2010)
  mm <- historicMeanAnnualBurnMap(polys, flammable, yrs)
  expect_equal(max(terra::values(mm), na.rm = TRUE), 1 / 20)
  ## unset: first to last year with a fire; unknown years (0) are not part of the record
  expect_identical(historicRecordYears(c(0L, 1995L, 2004L)), 1995:2004)
})

test_that("the fire regime counts each replicate, including one with no fire", {
  fs <- data.frame(rep = c(1, 1, 1, 2), year = c(5, 15, 15, 30), simSize = c(2, 10, 0, 4))
  r <- simulatedFireRegime(fs, reps = 1:3, simTimes = c(0, 40), pixelHa = 5)
  expect_identical(r$rep, 1:3)
  expect_identical(r$nFires, c(2L, 1L, 0L)) ## a fire of 0 cells is not a fire
  expect_equal(r$areaBurnedHaPerYear, c(12 * 5 / 40, 4 * 5 / 40, 0))
  expect_equal(r$largestFireHa, c(50, 20, 0))
})

test_that("the observed regime counts flammable land only", {
  flammable <- toyRaster(c(1, 1, NA, 1)) ## cell 3 cannot burn
  polys <- terra::vect(c("POLYGON ((0 0, 2 0, 2 1, 0 1, 0 0))", "POLYGON ((0 1, 1 1, 1 2, 0 2, 0 1))"),
                       crs = "EPSG:3857") ## the bottom row (cells 3, 4) and cell 1
  polys$YEAR <- c(2001L, 2002L)
  yrs <- 2001:2010
  historic <- terra::mask(historicMeanAnnualBurnMap(polys, flammable, yrs), flammable)
  o <- observedFireRegime(polys, historic, flammable, yrs)
  pixHa <- prod(terra::res(flammable)) / 1e4
  expect_equal(o$areaBurnedHaPerYear, 2 * pixHa / 10) ## cells 1 and 4, once each, over 10 years
  expect_equal(o$largestFireHa, pixHa)                ## each fire covers one flammable cell
  expect_identical(c(o$firstYear, o$lastYear), c(2001L, 2010L))
})

test_that("fire per block keeps blocks with no fire as zero", {
  fs <- data.frame(rep = c(1, 1, 2), year = c(5, 25, 6), simSize = c(1, 3, 2))
  b <- fireBlockTable(fs, reps = 1:2, simTimes = c(0, 30), blockYears = 10, pixelHa = 2)
  expect_equal(nrow(b), 2 * 3)
  expect_equal(b$nFires[b$rep == 1], c(1, 0, 1))
  expect_equal(b$areaBurnedHa[b$rep == 1], c(2, 0, 6))
  expect_equal(b$areaBurnedHa[b$rep == 2], c(4, 0, 0))
})
