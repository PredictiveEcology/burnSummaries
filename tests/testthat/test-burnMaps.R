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
