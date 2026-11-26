## updateTimeSinceFire(), the bookkeeping of the update_tsf event, on its own and against theory.

test_that("time since fire ages every cell by the timestep, resets burned cells and keeps NA", {
  expect_identical(updateTimeSinceFire(c(5L, NA, 0L, 7L), c(0, 0, 1, 1), 1L), c(6L, NA, 0L, 0L))
  expect_identical(updateTimeSinceFire(c(5L, NA, 0L, 7L), c(0, NA, 0, 1), 10L), c(15L, NA, 10L, 0L))
  ## a burn map read as a one-column matrix, as `raster[]` gives it
  expect_identical(updateTimeSinceFire(matrix(c(5L, 7L)), matrix(c(1, 0)), 1L), c(0L, 8L))
})

test_that("with the same chance of burning in every cell and year, ages follow the negative exponential", {
  ## Van Wagner (1978): when every cell has the same chance p of burning each year, independently of its
  ## age, time since fire settles to a negative exponential, with a mean age equal to the fire cycle
  ## 1/p (Johnson and Van Wagner 1985). In whole years, with a cell reading 0 in the year it burns, ages
  ## are geometric, the discrete form (Boychuk et al. 1997, p. 150): the share of cells at least t years
  ## old is (1 - p)^t and the mean age (1 - p) / p. Independent cells
  ## are the case in which this holds for one landscape at one time; with fires that spread, or burn
  ## rates that vary between years, it holds only on average over time (Boychuk et al. 1997).
  withr::local_seed(20261125)
  p <- 0.02
  n <- 20000L
  tsf <- integer(n)
  for (yr in seq_len(600)) { ## (1 - p)^600 < 1e-5: hardly a cell keeps its starting age
    tsf <- updateTimeSinceFire(tsf, stats::rbinom(n, 1, p), 1L)
  }
  expect_equal(mean(tsf), (1 - p) / p, tolerance = 0.03)
  for (t in c(10, 50, 100)) {
    expect_equal(mean(tsf >= t), (1 - p)^t, tolerance = 0.05)
  }
})
